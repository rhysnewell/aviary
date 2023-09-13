import pysam

samfile = pysam.AlignmentFile(snakemake.input[0], 'rb')


cutoff = 0.05
mapped_full = set()
overlap = set()
partial_map = set()
unmapped = set()
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
        length = read.infer_read_length()
        if clipped_start/length <= cutoff and clipped_end/length <= cutoff:
            mapped_full.add(read.query_name)
        elif clipped_start/length <= cutoff and read.reference_end > read.reference_length - 100:
            overlap.add(read.query_name)
        elif clipped_end/length <= cutoff and read.reference_start < 100:
            overlap.add(read.query_name)
        else:
            partial_map.add(read.query_name)
    else:
        unmapped.add(read.query_name)


if snakemake.params[0] == 'no_full':
    all_reads = unmapped.union(partial_map).union(overlap)
    out_set = all_reads - mapped_full
with open(snakemake.output[0], 'w') as o:
    for i in out_set:
        o.write(i + '\n')