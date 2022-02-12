import pysam
import os

contig_bins = {}
outlength = {}
samfile = pysam.AlignmentFile(snakemake.input.long_bam, 'rb')
gotten_contigs = set()
maxbin = 0
outreads = {}
outbases = {}
outreads500 = {}
outbases500 = {}

for bins in os.listdir(snakemake.input.metabat_done[:-4]):
    if not bins.startswith(
            'binned_contigs') or bins == "binned_contigs.unbinned" or bins == "binned_contigs.lowDepth" or bins == "binned_contigs.tooShort":
        continue
    bin = bins.split('.')[1]
    if int(bin) > maxbin:
        maxbin = int(bin)
    outlength[bin] = 0
    outreads[bin] = set()
    outbases[bin] = 0
    outreads500[bin] = set()
    outbases500[bin] = 0
    with open(os.path.join(snakemake.input.metabat_done[:-4], bins)) as f:
        for line in f:
            contig_bins[line.rstrip()] = bin
            gotten_contigs.add(line.rstrip())
            outlength[bin] += samfile.get_reference_length(line.rstrip())

cutoff = 0.05

for read in samfile.fetch(until_eof=True):
    start = True
    if not read.cigartuples is None:
        clipped_start = 0
        clipped_end = 0
        for i in read.cigartuples:
            if i[0] == 4 or i[0] == 5:
                if start:
                    clipped_start += i[1]
                else:
                    clipped_end += i[1]
            else:
                start = False
        if read.reference_name in contig_bins:
            bin = contig_bins[read.reference_name]
            length = read.infer_read_length()
            if (clipped_start / length <= cutoff and clipped_end / length <= cutoff) or \
                    (clipped_start / length <= cutoff and read.reference_end > read.reference_length - 100) or \
                    (clipped_end / length <= cutoff and read.reference_start < 100) or \
                    (read.reference_start < 100 and read.reference_end > read.reference_length - 100):
                outreads[bin].add(read.query_name)
                outbases[bin] += length

try:
    os.makedirs("data/binned_reads")
except FileExistsError:
    pass

out_dict = {}
for i in outreads:
    with open("data/binned_reads/r" + i + '.long.list', 'w') as read_list:
        reads = outreads[i]
        bases = outbases[i]
        for j in reads:
            read_list.write(j + '\n')
    out_dict[i] = [i, "data/binned_reads/r" + i + '.long.list', str(outlength[i]), str(bases)]

samfile = pysam.AlignmentFile(snakemake.input.short_bam, 'rb')
outreads = {}
outbases = {}
for read in samfile.fetch(until_eof=True):
    start = True
    if read.is_proper_pair:
        if read.reference_name in contig_bins:
            bin = contig_bins[read.reference_name]
            length = read.infer_read_length()
            if not bin in outreads:
                outreads[bin] = set()
                outbases[bin] = 0
            outreads[bin].add(read.query_name)
            outbases[bin] += length
for i in outreads:
    with open("data/binned_reads/r" + i + '.short.list', 'w') as read_list:
        for j in outreads[i]:
            if reads == 'data/short_reads.fastq.gz' or snakemake.config['short_reads_2'] == 'none':
                if j[:-2] in ['/1', '/2']:
                    read_list.write(j[:-2] + '/1\n')
                    read_list.write(j[:-2] + '/2\n')
                else:
                    read_list.write(j + '/1\n')
                    read_list.write(j + '/2\n')
            else:
                read_list.write(j + '\n')
    out_dict[i] += ["data/binned_reads/r" + i + '.short.list', str(outbases[i])]

with open(snakemake.output.list, 'w') as o:
    for i in sorted(list(out_dict), key=int):
        o.write('\t'.join(out_dict[i]) + '\n')
