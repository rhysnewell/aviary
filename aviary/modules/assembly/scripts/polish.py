import os
import sys
from subprocess import run, Popen, PIPE, STDOUT
import random
import shutil
import logging
import tempfile

def clean_short_reads(
    cat_or_zcat: str,
    read_path: str,
    read_pair: str,
    output_path: str,
    threads: int,
    log: str,
):
    cat_cmd = f"{cat_or_zcat} {read_path}".split()
    sed_cmd = f"""sed s/@/@{read_pair}_/""".split()

    with open(log, "a") as logf:
        logf.write(f"Shell command: {' '.join(cat_cmd)} | {' '.join(sed_cmd)} > {output_path}\n")
        logf.write(' '.join(sed_cmd))
        with open(output_path, 'a') as out:
            cat = Popen(cat_cmd, stdout=PIPE, stderr=logf)
            sed = Popen(sed_cmd, stdin=cat.stdout, stdout=out, stderr=logf)

            sed.wait()
            cat.wait()

            logf.write(f"cat return: {cat.returncode}\n")
            logf.write(f"sed return: {sed.returncode}\n")

def minimap2_process(
    minimap2_type: str,
    reference: str,
    reads: str,
    threads: int,
    output_paf: str,
    log: str,
):
    minimap2_cmd = f"minimap2 -x {minimap2_type} -t {threads} {reference} {reads}".split()

    with open(log, "a") as logf:
        with open(output_paf, 'w') as out:
            Popen(minimap2_cmd, stdout=out, stderr=logf).wait()

def run_seqkit(
    reads,
    pattern_file: str,
    output_file: str,
    threads: int,
    log: str,
):
    seqkit_cmd = f"seqkit -j {threads} grep  --pattern-file {pattern_file} {reads}".split()
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
            Popen(racon_cmd, stdout=out, stderr=logf).wait()


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
        # minimap2.stdout.close()
        # # minimap2.wait()
        # # samtools.stdin.close()
        # samtools.communicate()

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
    reference_filter: str,
    max_cov: int,
    illumina: bool,
    long_read_type: str,
    coassemble: bool,
    threads: int,
    log: str,
):
    # out = "data/polishing"

    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    random.seed(89)

    # Whether contigs are polished with illumina or long read
    if illumina:
        if reference_filter != 'none':
            reads = "data/short_reads.fastq.gz"
        elif short_reads_2 != 'none':

            # Racon can't handle paired end reads. It treats them as singled-ended. But when you have paired end reads
            # in separate files they can share the same read name, so we need to alter the read name based on the pair
            if not os.path.exists("data/short_reads.racon.1.fastq.gz"):
                for reads1, reads2 in zip(short_reads_1, short_reads_2):
                    cat_or_zcat1 = "cat"
                    cat_or_zcat2 = "cat"
                    if reads1[-3::] == ".gz":
                        cat_or_zcat1 = "zcat"
                    if reads2[-3::] == ".gz":
                        cat_or_zcat2 = "zcat"

                    clean_short_reads(cat_or_zcat1, reads1, 1, "data/short_reads.racon.1.fastq", threads, log)
                    clean_short_reads(cat_or_zcat2, reads2, 2, "data/short_reads.racon.1.fastq", threads, log)
                    if not coassemble:
                        break

                with open(log, "a") as logf:
                    run(f"pigz -p {threads} --fast data/short_reads.racon.1.fastq".split(), stdout=logf, stderr=STDOUT)

            pe1 = "data/short_reads.racon.1.fastq.gz"
            reads = [pe1]
        else:
            if len(short_reads_1) == 1 or not coassemble:
                pe1 = short_reads_1[0]
            else:
                if not os.path.exists("data/short_reads.1.fastq.gz"):
                    for reads1 in short_reads_1:
                        cat_or_zcat1 = "cat"
                        if reads1[-3::] == ".gz":
                            cat_or_zcat1 = "zcat"
                        
                        with open("data/short_reads.1.fastq", 'a') as out:
                            cat_cmd = f"{cat_or_zcat1} {reads1}".split()
                            with open(log, "a") as logf:
                                Popen(cat_cmd, stdout=out, stderr=logf).wait()

                    with open(log, "a") as logf:
                        run("pigz -p {threads} --fast data/short_reads.1.fastq".split(), stdout=logf, stderr=STDOUT)

                pe1 = "data/short_reads.1.fastq.gz"
            reads = [pe1]
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
                    if reads != "data/short_reads.fastq.gz":
                        minimap2_process("sr", reference, ' '.join(reads), threads, paf, log)
                    else:
                        minimap2_process("sr", reference, reads, threads, paf, log)
                elif long_read_type in ['ont', 'ont_hq']:
                    sys.exit("ONT reads are not supported for racon polishing")
                else:
                    minimap2_process("map-pb", reference, reads, threads, paf, log)

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
                    qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
                    # if illumina:
                    #     if qname[-2:] in ['/1', '/2']:
                    #         qname = qname[:-2]
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
                # if (reads == 'data/short_reads.fastq.gz' or short_reads_2 == 'none') and illumina:
                #     print("Using paired end reads and adding /1 and /2 to read names")
                # else:
                #     print("Using single end reads")

                for i in included_reads:
                    o.write(i + '\n')
                    # if (reads == 'data/short_reads.fastq.gz' or short_reads_2 == 'none') and illumina:
                    #     o.write(i + '/1\n')
                    #     o.write(i + '/2\n')
                    # else:
                    #     o.write(i + '\n')
            logging.info("Retrieving reads...")
            if not isinstance(reads, str):
                for read in reads:
                    pattern_file = f"{output_dir}/reads.{output_prefix}.{rounds}.lst"
                    output_file = f"{output_dir}/reads.{output_prefix}.{rounds}.fastq.gz"
                    run_seqkit(
                        read,
                        pattern_file,
                        output_file=output_file,
                        threads=threads,
                        log=log,
                    )

            else:
                pattern_file = f"{output_dir}/reads.{output_prefix}.{rounds}.lst"
                output_file = f"{output_dir}/reads.{output_prefix}.{rounds}.fastq.gz"
                run_seqkit(
                    reads,
                    pattern_file,
                    output_file=output_file,
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
        
        # bam = os.path.join(output_dir, 'alignment.%s.1.bam') % (output_prefix)
        # print("Generating BAM file: %s for medaka..." % (bam))
        # # we just run medaka once: https://twitter.com/rrwick/status/1158278701819125760
        # # Twitter is a valid source of information :) do not question.
        # run_minimap_with_samtools(
        #     reference,
        #     reads,
        #     output_file=bam,
        #     threads=threads,
        # )
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
        # os.remove("data/short_reads.racon.2.fastq.gz")


if __name__ == "__main__":

    short_reads_1 = snakemake.config["short_reads_1"]
    short_reads_2 = snakemake.config["short_reads_2"]
    input_fastq = snakemake.input.fastq
    reference = snakemake.input.fasta
    reference_filter = snakemake.config["reference_filter"]
    output_dir = "data/polishing"
    output_prefix = snakemake.params.prefix
    output_fasta = snakemake.output.fasta
    rounds = snakemake.params.rounds
    long_read_type = snakemake.config["long_read_type"]
    medaka_model = snakemake.config["medaka_model"]
    illumina = snakemake.params.illumina

    max_cov = snakemake.params.maxcov
    threads = snakemake.threads
    coassemble = snakemake.params.coassemble
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    run_polish(
        short_reads_1,
        short_reads_2,
        input_fastq,
        reference=reference,
        reference_filter=reference_filter,
        output_dir=output_dir,
        output_prefix=output_prefix,
        output_fasta=output_fasta,
        polishing_rounds=rounds,
        long_read_type=long_read_type,
        medaka_model=medaka_model,
        illumina=illumina,
        max_cov=max_cov,
        threads=threads,
        coassemble=coassemble,
        log=log,
    )