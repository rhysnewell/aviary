from subprocess import run, STDOUT
import os


def get_coverage(
    long_reads,
    short_reads_1,
    short_reads_2,
    long_read_type: str,
    input_fasta: str,
    tmpdir: str,
    threads: int,
    log: str,
):
    if tmpdir: os.environ["TMPDIR"] = tmpdir

    if long_reads != "none" and not os.path.exists("data/long_cov.tsv"):
        if long_read_type in ["ont", "ont_hq"]:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped --min-read-percent-identity 0.85 --output-file data/long_cov.tsv".split()
            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

        elif long_read_type in ["rs", "sq", "ccs", "hifi"]:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-pb -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped --min-read-percent-identity 0.9 --output-file data/long_cov.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

        else:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped --min-read-percent-identity 0.85 --output-file data/long_cov.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

    if short_reads_2 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -m metabat --bam-file-cache-directory data/binning_bams/ --discard-unmapped --output-file data/short_cov.tsv".split()
        with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

    elif short_reads_1  != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -m metabat --bam-file-cache-directory data/binning_bams/ --discard-unmapped --output-file data/short_cov.tsv".split()

        with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

    # Concatenate the two coverage files if both long and short exist
    if long_reads != "none" and (short_reads_1 != "none"):
        short_count = len(short_reads_1)
        long_count = len(long_reads)
        with open('data/coverm.cov', 'w') as file3:
            with open('data/short_cov.tsv', 'r') as file1:
                with open('data/long_cov.tsv', 'r') as file2:
                    for idx, (line1, line2) in enumerate(zip(file1, file2)):
                        long_values = line2.strip().split("\t")[1::]
                        del long_values[0::3] # delete length value
                        short_values = line1.strip().split("\t")
                        if idx != 0:
                            long_cov = sum([float(x) for x in long_values[0::2]])
                            short_cov = sum([float(x) for x in short_values[3::2]])
                            tot_depth = (long_cov + short_cov) / (short_count + long_count)
                            line = (
                                "{contig_name}\t"
                                "{contig_length}\t"
                                "{total_avg_depth}\t"
                                "{short}\t"
                                "{long}"
                            ).format(
                                contig_name = short_values[0],
                                contig_length = short_values[1],
                                total_avg_depth = tot_depth,
                                short = "\t".join(short_values[3::]),
                                long = "\t".join(long_values)
                            )
                            print(line, file=file3)
                        else:
                            line = (
                                "{short}\t{long}"
                            ).format(
                                short = "\t".join(short_values),
                                long = "\t".join(long_values)
                            )
                            print(line, file=file3)
                            
    elif long_reads != "none":
        with open('data/coverm.cov', 'w') as file3:
            with open('data/long_cov.tsv', 'r') as file1:
                for idx, line1 in enumerate(file1):
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
                            "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                        ).format(
                            samples = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
                            
    elif short_reads_1 != "none":  # rename shrot reads cov if only they exist
        os.rename("data/short_cov.tsv", "data/coverm.cov")

    if long_reads != "none":
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
                            "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                        ).format(
                            samples = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
            # os.remove("data/long_cov.tsv")

    try:
        os.makedirs("data/maxbin_cov/")
    except OSError:
        pass
    with open("data/coverm.cov") as f, open("data/maxbin.cov.list", "w") as o:
        cov_list = []
        contig_list = []
        f.readline()
        for line in f:
            contig = line.split()[0]
            contig_list.append(contig)
            depths = line.split()[3:]
            cov_list.append([])
            for i in range(0, len(depths), 2):
                cov_list[-1].append(depths[i])
        for i in range(len(cov_list[0])):
            with open("data/maxbin_cov/p%d.cov" % i, 'w') as oo:
                for j in range(len(contig_list)):
                    oo.write(contig_list[j] + '\t' + cov_list[j][i] + '\n')
            o.write("data/maxbin_cov/p%d.cov\n" % i)
    # with open("data/binning_bams/done", 'w') as o:
    #     o.write('done')


if __name__ == '__main__':
    long_reads = snakemake.config["long_reads"]
    short_reads_1 = snakemake.config["short_reads_1"]
    short_reads_2 = snakemake.config["short_reads_2"]
    long_read_type = snakemake.config["long_read_type"][0]
    input_fasta = snakemake.input.input_fasta
    tmpdir = snakemake.params.tmpdir
    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    get_coverage(
        long_reads,
        short_reads_1,
        short_reads_2,
        long_read_type,
        input_fasta,
        tmpdir,
        threads,
        log,
    )


