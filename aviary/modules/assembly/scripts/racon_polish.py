import os
import subprocess
import random
import shutil
import logging

out = "data/racon_polishing"
max_cov = snakemake.params.maxcov

try:
    os.makedirs(out)
except FileExistsError:
    pass

random.seed(89)
reference = snakemake.input.fasta

# Whether contigs are polished with illumina or long read
if snakemake.params.illumina:
    if snakemake.config['reference_filter'] != 'none':
        reads = "data/short_reads.fastq.gz"
    elif snakemake.config['short_reads_2'] != 'none':
        # if len(snakemake.config['short_reads_2']) == 1 or not snakemake.params.coassemble:
        #     pe1 = snakemake.config['short_reads_1'][0]
        #     pe2 = snakemake.config['short_reads_2'][0]
        # else:
        # Racon can't handle paired end reads. It treats them as singled-ended. But when you have paired end reads
        # in separate files they can share the same read name, so we need to alter the read name based on the pair
        if not os.path.exists("data/short_reads.racon.1.fastq.gz"):
            for reads1, reads2 in zip(snakemake.config['short_reads_1'], snakemake.config['short_reads_2']):
                cat_or_zcat1 = "cat"
                cat_or_zcat2 = "cat"
                if reads1[-3::] == ".gz":
                    cat_or_zcat1 = "zcat"
                if reads2[-3::] == ".gz":
                    cat_or_zcat2 = "zcat"

                subprocess.Popen(f"{cat_or_zcat1} {reads1} | sed 's/@/@1_/' >> data/short_reads.racon.1.fastq", shell=True).wait()
                subprocess.Popen(f"{cat_or_zcat2} {reads2} | sed 's/@/@2_/' >> data/short_reads.racon.1.fastq", shell=True).wait()
                if not snakemake.params.coassemble:
                    break

            subprocess.Popen(f"pigz -p {snakemake.threads} --fast data/short_reads.racon.1.fastq", shell=True).wait()
            # subprocess.Popen(f"pigz -p {snakemake.threads} --fast data/short_reads.racon.2.fastq", shell=True).wait()

        pe1 = "data/short_reads.racon.1.fastq.gz"
        # pe2 = "data/short_reads.racon.2.fastq.gz"
        reads = [pe1]
    else:
        if len(snakemake.config['short_reads_1']) == 1 or not snakemake.params.coassemble:
            pe1 = snakemake.config['short_reads_1'][0]
        else:
            if not os.path.exists("data/short_reads.1.fastq.gz"):
                for reads1 in snakemake.config['short_reads_1']:
                    cat_or_zcat1 = "cat"
                    if reads1[-3::] == ".gz":
                        cat_or_zcat1 = "zcat"
                    subprocess.Popen(f"{cat_or_zcat1} {reads1} >> data/short_reads.1.fastq", shell=True).wait()

                subprocess.Popen(f"pigz -p {snakemake.threads} --fast data/short_reads.1.fastq",
                                 shell=True).wait()
            pe1 = "data/short_reads.1.fastq.gz"
        reads = [pe1]
else:
    reads = snakemake.input.fastq

for rounds in range(snakemake.params.rounds):
    paf = os.path.join(out, 'alignment.%s.%d.paf') % (snakemake.params.prefix, rounds)
    print("Generating PAF file: %s for racon round %d..." % (paf, rounds))

    # Generate PAF mapping files
    if not os.path.exists(paf): # Check if mapping already exists
        if snakemake.params.illumina:
            if reads != "data/short_reads.fastq.gz":
                subprocess.Popen("minimap2 -t %d -x sr %s %s > %s" % (snakemake.threads, reference, ' '.join(reads), paf),
                                 shell=True).wait()
            else:
                subprocess.Popen("minimap2 -t %d -x sr %s %s > %s" % (snakemake.threads, reference, reads, paf),
                                 shell=True).wait()
        elif snakemake.config["long_read_type"] in ['ont', 'ont_hq']:
            subprocess.Popen("minimap2 -t %d -x map-ont %s %s > %s" % (snakemake.threads, reference, reads, paf), shell=True).wait()
        else:
            subprocess.Popen("minimap2 -t %d -x map-pb %s %s > %s" % (snakemake.threads, reference, reads, paf), shell=True).wait()

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
    with open(reference) as ref_file, open(os.path.join(out, "filtered.%s.%d.fa" % (snakemake.params.prefix, rounds)), 'w') as o:
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
    with open(paf) as f, open(os.path.join(out, "filtered.%s.%d.paf" % (snakemake.params.prefix, rounds)), 'w') as paf_file:
        for line in f:
            qname, qlen, qstart, qstop, strand, ref, rlen, rstart, rstop = line.split()[:9]
            qlen, qstart, qstop, rlen, rstart, rstop = map(int, [qlen, qstart, qstop, rlen, rstart, rstop])
            if snakemake.params.illumina:
                if qname[:-2] in ['/1', '/2']:
                    qname = qname[:-2]
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
    with open(os.path.join(out, "reads.%s.%d.lst" % (snakemake.params.prefix, rounds)), "w") as o:
        for i in included_reads:
            if (reads == 'data/short_reads.fastq.gz' or snakemake.config['short_reads_2'] == 'none') and snakemake.params.illumina:
                o.write(i + '/1\n')
                o.write(i + '/2\n')
            else:
                o.write(i + '\n')
    logging.info("Retrieving reads...")
    if not isinstance(reads, str):
        for read in reads:
            seqkit_command = f"seqkit -j {snakemake.threads} grep --pattern-file {out}/reads.{snakemake.params.prefix}.{rounds}.lst {read} | pigz -p {snakemake.threads} >> {out}/reads.{snakemake.params.prefix}.{rounds}.fastq.gz"
            print(seqkit_command)
            subprocess.Popen(seqkit_command, shell=True).wait()
    else:
        seqkit_command = f"seqkit -j {snakemake.threads} grep --pattern-file {out}/reads.{snakemake.params.prefix}.{rounds}.lst {reads} | pigz -p {snakemake.threads} >> {out}/reads.{snakemake.params.prefix}.{rounds}.fastq.gz"
        print(seqkit_command)
        subprocess.Popen(seqkit_command, shell=True).wait()

    print("Performing round %d of racon polishing..." % rounds)
    print("racon -m 8 -x -6 -g -8 -w 500 -t %d -u %s/reads.%s.%d.fastq.gz %s/filtered.%s.%d.paf %s/filtered.%s.%d.fa"
                     " > %s/filtered.%s.%d.pol.fa" % (snakemake.threads, out, snakemake.params.prefix, rounds, out,
                                                   snakemake.params.prefix, rounds, out, snakemake.params.prefix, rounds,
                                                   out, snakemake.params.prefix, rounds))
    subprocess.Popen("racon -m 8 -x -6 -g -8 -w 500 -t %d -u %s/reads.%s.%d.fastq.gz %s/filtered.%s.%d.paf %s/filtered.%s.%d.fa"
                     " > %s/filtered.%s.%d.pol.fa" % (snakemake.threads, out, snakemake.params.prefix, rounds, out,
                                                   snakemake.params.prefix, rounds, out, snakemake.params.prefix, rounds,
                                                   out, snakemake.params.prefix, rounds), shell=True).wait()

    with open(os.path.join(out, "combined.%s.%d.pol.fa" % (snakemake.params.prefix, rounds)), "w") as o:
        with open(os.path.join(out, "filtered.%s.%d.pol.fa" % (snakemake.params.prefix, rounds))) as f:
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
    reference = os.path.join(out, "combined.%s.%d.pol.fa" % (snakemake.params.prefix, rounds))


if os.path.exists("data/short_reads.racon.1.fastq.gz"):
    os.remove("data/short_reads.racon.1.fastq.gz")
    # os.remove("data/short_reads.racon.2.fastq.gz")
shutil.copy2(reference, snakemake.output.fasta)
