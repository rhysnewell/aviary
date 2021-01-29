import os

with open('data/long.cov', 'w') as file3:
    with open('data/long_cov.tsv', 'r') as file:
        for idx, line1 in enumerate(file):
            long_values = line1.strip().split("\t")
            del long_values[4::3] # delete extra length values
            long_count = len(long_values[2::2])
            if idx != 0:
                long_cov = sum([float(x) for x in long_values[2::2]])
                tot_depth = long_cov / long_count
                line = (
                    "{contig_name}\t"
                    "{contig_length}\t"
                    "{total_avg_depth}\t"
                    "{long}"
                ).format(
                    contig_name = long_values[0],
                    contig_length = long_values[1],
                    total_avg_depth = tot_depth,
                    long = "\t".join(long_values[2::])
                )
                print(line, file=file3)
            else:
                line = (
                    "{contig_name}\t{contig_length}\ttotalAvgDepth\t{samples}"
                ).format(
                    contig_name = 'contigName',
                    contig_length = 'contigLen',
                    samples = "\t".join(long_values[2::])
                )
                print(line, file=file3)
    # os.remove("data/long_cov.tsv")

