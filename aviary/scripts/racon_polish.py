import sys
import os
import subprocess
import random
import shutil
import gzip

out = "data/racon_polishing"
max_cov = snakemake.params.maxcov

try:
    os.makedirs(out)
except FileExistsError:
    pass

random.seed(89)
reference = snakemake.input.fasta

if snakemake.params.illumina:
    # with gzip.open(reads, 'rt') as f:
    #     line1 = f.readline()
    #     f.readline()
    #     f.readline()
    #     f.readline()
    #     line5 = f.readline()
    # if line1.split()[0] == line5.split()[0]:
    #     subprocess.Popen('zcat ' + reads + ' | awk \'{{if (NR % 8 == 1) {{print $1 "/1"}} else if (NR % 8 == 5) {{print $1 "/2"}} ' \
    #                      'else if (NR % 4 == 3){{print "+"}} else {{print $0}}}}\' | gzip > data/short_reads_racon.fastq.gz', shell=True).wait()
    if snakemake.config['reference_filter'] != 'none':
        reads = "data/short_reads.fastq.gz"
    elif snakemake.config['short_reads_2'] != 'none':
        reads = [' '.join([pe1, pe2]) for pe1, pe2 in zip(snakemake.config['short_reads_1'], snakemake.config['short_reads_2'])]
    else:
        reads = snakemake.config['short_reads_1']
else:
    reads = snakemake.input.fastq

for rounds in range(snakemake.params.rounds):
    paf = os.path.join(out, 'alignment.%s.%d.paf') % (snakemake.params.prefix, rounds)
    if snakemake.params.illumina:
        if reads != "data/short_reads.fastq.gz":
            subprocess.Popen("minimap2 -t %d -x sr %s %s > %s" % (snakemake.threads, reference, ' '.join(reads), paf),
                             shell=True).wait()
        else:
            subprocess.Popen("minimap2 -t %d -x sr %s %s > %s" % (snakemake.threads, reference, reads, paf),
                             shell=True).wait()
    elif config["long_read_type"] == 'ont':
        subprocess.Popen("minimap2 -t %d -x map-ont %s %s > %s" % (snakemake.threads, reference, reads, paf), shell=True).wait()
    else:
        subprocess.Popen("minimap2 -t %d -x map-pb %s %s > %s" % (snakemake.threads, reference, reads, paf), shell=True).wait()

    cov_dict = {}
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
                qname = qname[:-2]
            if ref in low_cov:
                paf_file.write(line)
                included_reads.add(qname)
            elif ref in high_cov:
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
            # if reads == 'data/short_reads.fastq.gz':
            #     o.write(i + '/1\n')
            #     o.write(i + '/2\n')
            # else:
            o.write(i + '\n')
    subprocess.Popen("seqkit grep --pattern-file %s/reads.%s.%d.lst | gzip > %s/reads.%s.%d.fastq.gz" % \
                     (out, snakemake.params.prefix, rounds, reads, out, snakemake.params.prefix, rounds), shell=True).wait()
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

shutil.copy2(reference, snakemake.output.fasta)