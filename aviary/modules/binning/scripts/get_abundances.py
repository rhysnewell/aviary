from subprocess import run, STDOUT
import os

def run_coverm(
    reads: str,
    minimap2_type: str,
    output_file: str,
    read_type: str,
    threads: int,
    strain_analysis: bool,
    output_dir: str,
    log: str,
):
    strain_analysis_flag = f"--bam-file-cache-directory {output_dir} --discard-unmapped" if strain_analysis else ""

    coverm_cmd = f"coverm genome -t {threads} {strain_analysis_flag} -d bins/final_bins/ -m relative_abundance covered_fraction {read_type} {reads} -p {minimap2_type} --output-file {output_file} --min-covered-fraction 0.0 -x fna".split()

    with open(log, "a") as logf:
        run(coverm_cmd, stdout=logf, stderr=STDOUT)

def get_abundances(
    long_reads,
    short_reads_1,
    short_reads_2,
    long_read_type: str,
    threads: int,
    strain_analysis: bool,
    log: str,
):
    if long_reads != "none":
        if long_read_type in ["ont", "ont_hq"]:
            run_coverm(
                reads=" ".join(long_reads),
                minimap2_type="minimap2-ont",
                output_file="data/long_abundances.tsv",
                read_type="--single",
                threads=threads,
                strain_analysis=strain_analysis,
                output_dir="data/reads_mapped_to_mags/long/",
                log=log,
            )

        elif long_read_type in ["rs", "sq", "ccs", "hifi"]:
            run_coverm(
                reads=" ".join(long_reads),
                minimap2_type="minimap2-pb",
                output_file="data/long_abundances.tsv",
                read_type="--single",
                threads=threads,
                strain_analysis=strain_analysis,
                output_dir="data/reads_mapped_to_mags/long/",
                log=log,
            )

        else:
            run_coverm(
                reads=" ".join(long_reads),
                minimap2_type="minimap2-ont",
                output_file="data/long_abundances.tsv",
                read_type="--single",
                threads=threads,
                strain_analysis=strain_analysis,
                output_dir="data/reads_mapped_to_mags/long/",
                log=log,
            )

    if short_reads_2 != 'none':
        reads_str = []
        for r1, r2 in zip(short_reads_1, short_reads_2):
            reads_str.append(f"{r1} {r2}")
        reads_str = " ".join(reads_str)

        run_coverm(
            reads=reads_str,
            minimap2_type="minimap2-sr",
            output_file="data/short_abundances.tsv",
            read_type="--coupled",
            threads=threads,
            strain_analysis=strain_analysis,
            output_dir="data/reads_mapped_to_mags/short/",
            log=log,
        )

    elif short_reads_1 != 'none':
        run_coverm(
            reads=" ".join(short_reads_1),
            minimap2_type="minimap2-sr",
            output_file="data/short_abundances.tsv",
            read_type="--interleaved",
            threads=threads,
            strain_analysis=strain_analysis,
            output_dir="data/reads_mapped_to_mags/short/",
            log=log,
        )

    # Concatenate the two coverage files if both long and short exist
    if long_reads != "none" and short_reads_1 != "none":
        with open('data/coverm_abundances.tsv', 'w') as file3:
            with open('data/short_abundances.tsv', 'r') as file1:
                with open('data/long_abundances.tsv', 'r') as file2:
                    for line1, line2 in zip(file1, file2):
                        long_cov_line = "\t".join([l.strip() for l in line2.strip().split('\t')[1::]])
                        print(line1.strip(), "\t", long_cov_line, file=file3)
    elif long_reads != "none":  # rename long reads cov if only it exists
        os.rename("data/long_abundances.tsv", "data/coverm_abundances.tsv")
    elif short_reads_1 != "none":  # rename shrot reads cov if only they exist
        os.rename("data/short_abundances.tsv", "data/coverm_abundances.tsv")


if __name__ == '__main__':
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    get_abundances(
        snakemake.config["long_reads"],
        snakemake.config["short_reads_1"],
        snakemake.config["short_reads_2"],
        snakemake.config["long_read_type"][0],
        snakemake.threads,
        snakemake.config["strain_analysis"],
        log,
    )
