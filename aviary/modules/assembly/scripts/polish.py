#!/usr/bin/env python3
import os
import sys
import argparse
from subprocess import run, Popen, PIPE, STDOUT
import random
import shutil
import logging
import tempfile

# Snakemake runs this script standalone inside the `polishing` pixi env, where
# aviary is not pip-installed, so the package has to be put on sys.path by path
# rather than imported as an installed dependency. Matches get_coverage.py /
# get_abundances.py / fraction_recovered.py.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
from aviary import aligner_preset_args

# Rewrites each fastq header to "@<mate>_<id>", where id is the first header
# token minus any /1 or /2. mate=0 numbers an interleaved file's records 1, 2,
# 1, 2, ... (on single-end input that just alternates prefixes, harmlessly).
_RACON_READ_NAME_AWK = (
    'NR % 4 == 1 { id = substr($1, 2); sub(/\\/[12]$/, "", id); '
    'print "@" (mate ? mate : (NR % 8 == 1 ? 1 : 2)) "_" id; next } { print }'
)


def prepare_short_reads(reads_and_mates, output_path: str, threads: int, log: str) -> None:
    """Write short reads for racon with one plain, unique ID per read.

    racon treats every read as single-end and looks reads up by ID, so mates
    sharing an ID make it pair one mate's overlap with the other's sequence
    (isovic/racon#68): "unequal lengths in sequence and overlap file", or
    silently wrong polishing when the mates are the same length. Giving every
    read a suffix-free ID up front also means every aligner's PAF query names
    match these IDs exactly, so nothing downstream needs to reconcile them.

    reads_and_mates: (fastq path, mate) pairs -- mate 1 or 2 for an R1/R2
    file, 0 for an interleaved or single-end file.
    """
    if os.path.exists(output_path):
        return
    tmp_path = output_path + ".tmp"
    with open(log, "a") as logf, open(tmp_path, "wb") as out:
        pigz = Popen(["pigz", "-p", str(threads), "--fast"], stdin=PIPE, stdout=out, stderr=logf)
        for path, mate in reads_and_mates:
            logf.write(f"Preparing {path} (mate {mate}) for racon -> {output_path}\n")
            logf.flush()
            cat = Popen(["zcat" if path.endswith(".gz") else "cat", path], stdout=PIPE, stderr=logf)
            awk = Popen(["awk", "-v", f"mate={mate}", _RACON_READ_NAME_AWK],
                        stdin=cat.stdout, stdout=pigz.stdin, stderr=logf)
            cat.stdout.close()
            if awk.wait() != 0 or cat.wait() != 0:
                raise RuntimeError(f"preparing {path} for racon failed; see {log} for details")
        pigz.stdin.close()
        if pigz.wait() != 0:
            raise RuntimeError(f"compressing {output_path} failed; see {log} for details")
    os.replace(tmp_path, output_path)


# PAF-emitting command per short-read aligner. racon needs PAF, and each
# aligner spells that differently: minimap2/rammap take a preset with -x,
# strobealign's -x is a boolean switch, and minibwa needs the "map"
# subcommand with -f. Keys are the CoverM -p values that --short-read-mapper
# resolves to (including a --short-read-mapper-model suffix, if one was
# given), so this stays in step with the coverage side. bwa-mem/bwa-mem2 have
# no PAF mode at all and are handled separately by bwa_mem_paf_cmd() below.
def short_read_paf_cmd(mapper: str, reference: str, reads: str, threads: int) -> str:
    if mapper == "strobealign":
        return f"strobealign -x -t {threads} {reference} {reads}"
    for family in ("minimap2", "rammap"):
        if mapper.startswith(f"{family}-"):
            # The model half of the CoverM name is not the aligner's own -x
            # value (e.g. CoverM "lr-hq" is minimap2's "lr:hq"), so it has to be
            # translated rather than passed through -- see aligner_preset_args.
            preset = " ".join(aligner_preset_args(mapper.removeprefix(f"{family}-")))
            preset = f"{preset} " if preset else ""
            return f"{family} {preset}-t {threads} {reference} {reads}"
    if mapper == "minibwa":
        return f"minibwa map -f -x sr -t {threads} {reference} {reads}"
    raise ValueError(f"No PAF mode known for short-read mapper {mapper!r}")


# minibwa's PAF mode for long reads, mirroring short_read_paf_cmd()'s minibwa
# branch above but with -x lr instead of -x sr. -x belongs to `map`'s own
# option loop, not a top-level flag -- calling it minimap2-style as
# `minibwa -x lr ...` fails with "unknown command '-x'". Verified against the
# installed minibwa 0.6 binary end-to-end: `minibwa map -f -x lr` against
# test/data/pbsim.fq.gz mapped all 1869 reads to test/data/assembly.fasta,
# and the resulting PAF fed straight into racon without modification.
def minibwa_long_read_paf_cmd(reference: str, reads: str, threads: int) -> str:
    return f"minibwa map -f -x lr -t {threads} {reference} {reads}"


def index_reference(mapper: str, reference: str, log: str) -> None:
    """bwa/bwa-mem2/minibwa require a prebuilt on-disk index before `mem`/`map`
    can run, unlike strobealign/minimap2/rammap which index inline from the
    FASTA path in one invocation. Verified against the installed minibwa 0.6
    binary: `minibwa map -f -x sr <fasta> <reads>` aborts with "failed to load
    the index. ABORT!" because it treats its first positional argument as an
    index path, not a FASTA -- `minibwa index <fasta>` must run first, writing
    `<fasta>.l2b`/`<fasta>.mbw` alongside it, after which the original FASTA
    path doubles as the index argument to `map` (same convention as bwa). The
    reference is a fresh FASTA every racon round (see run_polish()), so this
    runs once per round.
    """
    index_bin = {"bwa-mem": "bwa", "bwa-mem2": "bwa-mem2", "minibwa": "minibwa"}[mapper]
    index_cmd = [index_bin, "index", reference]
    with open(log, "a") as logf:
        print("Running command:", " ".join(index_cmd), file=logf)
        returncode = run(index_cmd, stdout=logf, stderr=STDOUT).returncode
        logf.write(f"{index_bin} index return: {returncode}\n")
        if returncode != 0:
            raise RuntimeError(
                f"{index_bin} index failed (exit {returncode}) on {reference}; see {log} for details"
            )


def bwa_mem_paf_cmd(
    mapper: str, reference: str, reads: str, threads: int, output_paf: str, log: str,
) -> None:
    """bwa mem / bwa-mem2 mem only emit SAM, unlike every other supported
    short-read aligner which can be asked for PAF directly. Convert with
    paftools.js sam2paf, which ships as part of the polishing env's minimap2
    package, so no extra dependency is needed for the conversion step.
    """
    index_reference(mapper, reference, log)
    mem_bin = "bwa" if mapper == "bwa-mem" else "bwa-mem2"
    mem_cmd = [mem_bin, "mem", "-t", str(threads), reference, *reads.split()]
    sam2paf_cmd = ["paftools.js", "sam2paf", "-P", "-"]
    with open(log, "a") as logf:
        print("Running command:", " ".join(mem_cmd), "|", " ".join(sam2paf_cmd), file=logf)
        with open(output_paf, "w") as out:
            mem = Popen(mem_cmd, stdout=PIPE, stderr=logf)
            sam2paf = Popen(sam2paf_cmd, stdin=mem.stdout, stdout=out, stderr=logf)
            sam2paf.wait()
            mem.wait()
        logf.write(f"{mem_bin} return: {mem.returncode}\n")
        logf.write(f"sam2paf return: {sam2paf.returncode}\n")
        if mem.returncode != 0 or sam2paf.returncode != 0:
            raise RuntimeError(
                f"{mem_bin}/sam2paf pipeline failed (mem exit {mem.returncode}, "
                f"sam2paf exit {sam2paf.returncode}) generating {output_paf}; see {log} for details"
            )


def minimap2_process(
    minimap2_type: str,
    reference: str,
    reads: str,
    threads: int,
    output_paf: str,
    log: str,
    mapper: str = "minimap2",
):
    # rammap is a minimap2-compatible Rust implementation: same -x presets,
    # same PAF output, which this script parses directly further down. Only
    # the binary name differs, so --long-read-mapper selects it here too.
    #
    # rammap's default (non-CIGAR) PAF path derives target-start from the raw
    # seed anchor without clamping to 0, so an alignment landing at the very
    # start of the target underflows to ~u64::MAX instead of 0 and segfaults
    # racon downstream. minimap2 never hits this because its equivalent
    # anchor-extrapolation step (mm_align1 in align.c) always clamps negative
    # starts to 0, unconditionally. -c forces rammap through its full
    # base-level alignment path, which computes the same clamped coordinate
    # minimap2 always does. Confirmed via minimal repro: with -c, the two
    # previously-wrapped coordinates come back as 0/1 and racon completes.
    extra_args = ["-c"] if mapper == "rammap" else []
    minimap2_cmd = [mapper, "-x", minimap2_type, "-t", str(threads), *extra_args, reference, reads]

    with open(log, "a") as logf:
        # Logged so the aligner actually used is recoverable from the run.
        print("Running command:", " ".join(minimap2_cmd), file=logf)
        with open(output_paf, 'w') as out:
            returncode = Popen(minimap2_cmd, stdout=out, stderr=logf).wait()
        logf.write(f"{mapper} return: {returncode}\n")
        if returncode != 0:
            raise RuntimeError(
                f"{mapper} failed (exit {returncode}) generating {output_paf}; "
                f"see {log} for details"
            )

def run_seqkit(
    reads,
    pattern_file: str,
    output_file: str,
    threads: int,
    log: str,
):
    seqkit_cmd = ["seqkit", "-j", str(threads), "grep", "--pattern-file", pattern_file, reads]
    pigz_cmd = f"pigz -p {threads}".split()

    with open(log, "a") as logf:
        logf.write(f"Shell style: {' '.join(seqkit_cmd)} | {' '.join(pigz_cmd)} > {output_file}\n")

        with open(output_file, 'a') as out:
            seqkit = Popen(seqkit_cmd, stdout=PIPE, stderr=logf)
            pigz = Popen(pigz_cmd, stdin=seqkit.stdout, stdout=out, stderr=logf)
            pigz.wait()
            seqkit.wait()
            logf.write(f"seqkit return: {seqkit.returncode}\n")
            logf.write(f"pigz return: {pigz.returncode}\n")

def run_racon(
    reads: str,
    paf: str,
    reference: str,
    output_file: str,
    threads: int,
    log: str,
):
    racon_cmd = f"racon -m 8 -x -6 -g -8 -w 500 -t {threads} -u {reads} {paf} {reference}".split()

    with open(log, "a") as logf:
        logf.write(' '.join(racon_cmd))
        with open(output_file, 'w') as out:
            returncode = Popen(racon_cmd, stdout=out, stderr=logf).wait()
        logf.write(f"\nracon return: {returncode}\n")
        if returncode != 0:
            raise RuntimeError(
                f"racon failed (exit {returncode}) on {reads} against {paf}/{reference}; "
                f"see {log} for details"
            )


def run_minimap_with_samtools(
    reference: str,
    reads: str,
    threads: int,
    output_file: str,
    log: str,
):

    # write minimap2 output to temporary file

    with open(log, "a") as logf:
        minimap2_cmd = f"minimap2 -ax map-ont -t {threads} {reference} {reads}".split()
        logf.write(' '.join(minimap2_cmd))
        samtools_cmd = f"samtools view -F 4 -b -@ {threads-1} -o {output_file}".split()
        logf.write(' '.join(samtools_cmd))
        minimap2 = Popen(minimap2_cmd, stdout=PIPE, stderr=logf)
        samtools = Popen(samtools_cmd, stdin=minimap2.stdout, stderr=logf)
        samtools.wait()
        minimap2.wait()

    # check if output file exists and is not empty
    with open(log, "a") as logf:
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logf.write(f"{output_file} created.\n")
            if samtools.returncode == 0:
                logf.write("samtools successfully created bam file.\n")
            else:
                logf.write("samtools failed to create bam file.\n")
                logf.write(f"samtools return: {samtools.returncode}\n")
            return True
        else:
            logf.write(f"Error: {output_file} is empty or does not exist.\n")
            logf.write(f"samtools return: {samtools.returncode}\n")
            return False


def run_polish(
    short_reads_1,
    short_reads_2,
    input_fastq,
    output_dir: str,
    output_prefix: str,
    output_fasta: str,
    polishing_rounds: int,
    medaka_model: str,
    reference: str,
    max_cov: int,
    illumina: bool,
    long_read_type: str,
    coassemble: bool,
    threads: int,
    log: str,
    long_read_mapper: str = "rammap",
    short_read_mapper: str = "strobealign",
):
    if not reference or not os.path.exists(reference) or os.path.getsize(reference) == 0:
        os.makedirs(os.path.dirname(output_fasta) or ".", exist_ok=True)
        with open(log, "a") as logf:
            logf.write("Reference assembly is empty; skipping polishing.\n")
        if reference and os.path.exists(reference):
            shutil.copyfile(reference, output_fasta)
        else:
            with open(output_fasta, "w"):
                pass
        return

    # out = "data/polishing"

    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    random.seed(89)

    # Whether contigs are polished with illumina or long read
    if illumina:
        # Separate R1/R2 files, or one interleaved/single-end file per sample
        # (short_reads_2 == 'none', e.g. aviary's own QC output). Mates need
        # unique IDs either way -- see prepare_short_reads().
        if short_reads_2 != 'none':
            samples = [((r1, 1), (r2, 2)) for r1, r2 in zip(short_reads_1, short_reads_2)]
        else:
            samples = [((r1, 0),) for r1 in short_reads_1]
        if not coassemble:
            samples = samples[:1]
        reads = "data/short_reads.racon.1.fastq.gz"
        prepare_short_reads([f for sample in samples for f in sample], reads, threads, log)
    else:
        reads = input_fastq

    # use racon when using illumina or pacbio data
    if illumina or long_read_type not in ['ont', 'ont_hq']:
        for rounds in range(polishing_rounds):
            paf = os.path.join(output_dir, 'alignment.%s.%d.paf') % (output_prefix, rounds)
            with open(log, "a") as logf:
                logf.write("Generating PAF file: %s for racon round %d...\n" % (paf, rounds))

            # Generate PAF mapping files
            if not os.path.exists(paf): # Check if mapping already exists
                if illumina:
                    if short_read_mapper in ("bwa-mem", "bwa-mem2"):
                        bwa_mem_paf_cmd(
                            short_read_mapper, reference, reads, threads, paf, log)
                    else:
                        if short_read_mapper == "minibwa":
                            index_reference(short_read_mapper, reference, log)
                        cmd = short_read_paf_cmd(
                            short_read_mapper, reference, reads, threads)
                        with open(log, "a") as logf:
                            logf.write(f"Short-read PAF command: {cmd}\n")
                            logf.flush()
                            with open(paf, 'w') as out:
                                returncode = Popen(cmd.split(), stdout=out, stderr=logf).wait()
                            logf.write(f"{short_read_mapper} return: {returncode}\n")
                            if returncode != 0:
                                raise RuntimeError(
                                    f"{short_read_mapper} failed (exit {returncode}) generating "
                                    f"{paf}; see {log} for details"
                                )
                elif long_read_type in ['ont', 'ont_hq']:
                    sys.exit("ONT reads are not supported for racon polishing")
                elif long_read_mapper == "minibwa":
                    # minibwa takes a `map` subcommand rather than minimap2's
                    # bare `-x <preset>`, so minimap2_process() below cannot be
                    # reused for it (that would call `minibwa -x map-pb ...`
                    # and get "unknown command '-x'"). See
                    # minibwa_long_read_paf_cmd() for the working invocation.
                    index_reference(long_read_mapper, reference, log)
                    cmd = minibwa_long_read_paf_cmd(reference, reads, threads)
                    with open(log, "a") as logf:
                        logf.write(f"Long-read PAF command: {cmd}\n")
                        logf.flush()
                        with open(paf, 'w') as out:
                            returncode = Popen(cmd.split(), stdout=out, stderr=logf).wait()
                        logf.write(f"{long_read_mapper} return: {returncode}\n")
                        if returncode != 0:
                            raise RuntimeError(
                                f"{long_read_mapper} failed (exit {returncode}) generating "
                                f"{paf}; see {log} for details"
                            )
                else:
                    minimap2_process("map-pb", reference, reads, threads, paf, log,
                                     mapper=long_read_mapper)

            # minibwa's -f PAF output includes an explicit row per unmapped
            # read (all numeric/ref fields "*"), unlike strobealign/minimap2/
            # rammap/bwa-mem/bwa-mem2, which omit unmapped reads from PAF
            # output entirely. Strip those rows here, right after the PAF is
            # written, so every reader of this file -- both below in this
            # function and get_high_cov_contigs in assembly.smk, which reads
            # this same on-disk file directly -- can assume every line is a
            # real alignment without needing its own "*" guard. A no-op for
            # every other mapper, since they never emit these rows.
            filtered_lines = []
            with open(paf) as f:
                for line in f:
                    fields = line.split()[:9]
                    if len(fields) > 5 and fields[5] == '*':
                        continue
                    filtered_lines.append(line)
            tmp_paf = paf + '.tmp'
            with open(tmp_paf, 'w') as f:
                f.writelines(filtered_lines)
            os.replace(tmp_paf, paf)

            cov_dict = {}
            # Populate coverage dictionary,
            with open(paf) as f:
                for line in f:
                    qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
                    qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
                    if ref in cov_dict:
                        cov_dict[ref] += (rstop - rstart) / rlen
                    else:
                        cov_dict[ref] = (rstop - rstart) / rlen

            high_cov = set()
            low_cov = set()
            for i in cov_dict:
                if cov_dict[i] >= max_cov:
                    high_cov.add(i)
                else:
                    low_cov.add(i)

            no_cov = set()
            with open(reference) as ref_file, open(os.path.join(output_dir, "filtered.%s.%d.fa" % (output_prefix, rounds)), 'w') as o:
                for line in ref_file:
                    if line.startswith('>'):
                        name = line.split()[0][1:]
                        if name in low_cov or name in high_cov:
                            o.write(line)
                            getseq = True
                        else:
                            no_cov.add(name)
                            getseq = False
                    elif getseq:
                        o.write(line)

            included_reads = set()
            excluded_reads = set()
            with open(paf) as f, open(os.path.join(output_dir, "filtered.%s.%d.paf" % (output_prefix, rounds)), 'w') as paf_file:
                for line in f:
                    qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
                    if ref in low_cov:
                        paf_file.write(line)
                        included_reads.add(qname)
                    elif ref in high_cov:
                        # Down sample reads from high coverage contigs
                        sample_rate = max_cov / cov_dict[ref]
                        if qname in excluded_reads:
                            pass
                        elif qname in included_reads:
                            paf_file.write(line)
                        elif random.random() < sample_rate:
                            included_reads.add(qname)
                            paf_file.write(line)
                        else:
                            excluded_reads.add(qname)
            with open(os.path.join(output_dir, "reads.%s.%d.lst" % (output_prefix, rounds)), "w") as o:
                for i in included_reads:
                    o.write(i + '\n')
            logging.info("Retrieving reads...")
            run_seqkit(
                reads,
                f"{output_dir}/reads.{output_prefix}.{rounds}.lst",
                output_file=f"{output_dir}/reads.{output_prefix}.{rounds}.fastq.gz",
                threads=threads,
                log=log,
            )

            with open(log, "a") as logf:
                logf.write("Performing round %d of racon polishing..." % rounds)

            reads = f"{output_dir}/reads.{output_prefix}.{rounds}.fastq.gz"
            paf_file = f"{output_dir}/filtered.{output_prefix}.{rounds}.paf"
            reference = f"{output_dir}/filtered.{output_prefix}.{rounds}.fa"
            output_file = f"{output_dir}/filtered.{output_prefix}.{rounds}.pol.fa"

            run_racon(
                reads,
                paf_file,
                reference,
                output_file,
                threads=threads,
                log=log,
            )

            with open(os.path.join(output_dir, "combined.%s.%d.pol.fa" % (output_prefix, rounds)), "w") as o:
                with open(os.path.join(output_dir, "filtered.%s.%d.pol.fa" % (output_prefix, rounds))) as f:
                    gotten_set = set()
                    for line in f:
                        if line.startswith('>'):
                            gotten_set.add(line.split()[0][1:])
                        o.write(line)
                with open(reference) as f:
                    for line in f:
                        if line.startswith('>'):
                            name = line.split()[0][1:]
                            if name in gotten_set:
                                get_line = False
                            else:
                                get_line = True
                        if get_line:
                            o.write(line)
            reference = os.path.join(output_dir, "combined.%s.%d.pol.fa" % (output_prefix, rounds))

        # Copy the final polished reference to the output directory
        shutil.copyfile(reference, output_fasta)
    else:
        # polishing will be done by medaka
        if long_read_type not in ['ont', 'ont_hq']:
            sys.exit("ERROR: long_read_type must be ont or ont_hq for medaka polishing")
        
        with open(log, "a") as logf:
            logf.write("Running medaka...\n")
            medaka_cmd = f"medaka_consensus -t {threads} -i {reads} -m {medaka_model} -o data/polishing/ -d {reference}".split()
            logf.write(' '.join(medaka_cmd))
            run(medaka_cmd, stdout=logf, stderr=STDOUT)

        # copy the output to the expected location
        shutil.copyfile("data/polishing/consensus.fasta", output_fasta)
        # remove the intermediate files
        os.remove("data/polishing/calls_to_draft.bam")
        os.remove("data/polishing/calls_to_draft.bam.bai")
        os.remove("data/polishing/consensus_probs.hdf")


    if os.path.exists("data/short_reads.racon.1.fastq.gz"):
        os.remove("data/short_reads.racon.1.fastq.gz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Polish reads using racon or medaka.')
    parser.add_argument('--short-reads-1', nargs='*', help='Short reads 1')
    parser.add_argument('--short-reads-2', nargs='+', default='none', help='Short reads 2')
    parser.add_argument('--input-fastq', help='Input fastq file')
    parser.add_argument('--reference', help='Reference fasta file')
    parser.add_argument('--output-dir', default='data/polishing', help='Output directory')
    parser.add_argument('--output-prefix', help='Output prefix')
    parser.add_argument('--output-fasta', help='Output fasta file')
    parser.add_argument('--rounds', type=int, help='Number of polishing rounds')
    parser.add_argument('--long-read-type', help='Long read type')
    parser.add_argument('--long-read-mapper', default='rammap',
                        help='Aligner for long reads (rammap or minimap2)')
    parser.add_argument('--short-read-mapper', default='strobealign',
                        help='CoverM -p value for short reads; also selects the '
                             'aligner used for racon PAF generation')
    parser.add_argument('--medaka-model', help='Medaka model')
    parser.add_argument('--illumina', type=lambda x: x.lower() == 'true', nargs='?', const=True, default=False, help='Use illumina reads')
    parser.add_argument('--max-cov', type=int, default=100, help='Maximum coverage')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--coassemble', type=lambda x: x.lower() == 'true', nargs='?', const=True, default=False, help='Co-assemble')
    parser.add_argument('--log', default='polish.log', help='Log file')
    
    args = parser.parse_args()
    
    with open(args.log, "w") as logf:
        pass

    read2 = args.short_reads_2
    if read2 == ['none']:
        read2 = 'none'
    
    run_polish(
        args.short_reads_1,
        read2,
        args.input_fastq,
        reference=args.reference,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        output_fasta=args.output_fasta,
        polishing_rounds=args.rounds,
        long_read_type=args.long_read_type,
        long_read_mapper=args.long_read_mapper,
        short_read_mapper=args.short_read_mapper,
        medaka_model=args.medaka_model,
        illumina=args.illumina,
        max_cov=args.max_cov,
        threads=args.threads,
        coassemble=args.coassemble,
        log=args.log,
    )
