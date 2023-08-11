from subprocess import run
import os



def get_fraction_recovered(
    long_reads,
    short_reads_1,
    short_reads_2,
    input_fasta: str,
    long_read_type: str,
    threads: int,
):
    if long_reads != "none" and not os.path.exists("data/long_cov.tsv"):
        if long_read_type in ["ont", "ont_hq"]:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont --min-read-percent-identity 0.85 -o www/fraction_recovered/long_fraction_recovered.tsv".split()
            run(coverm_cmd)

        elif long_read_type in ["rs", "sq", "ccs", "hifi"]:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-pb --min-read-percent-identity 0.9 -o www/fraction_recovered/long_fraction_recovered.tsv".split()

            run(coverm_cmd)

        else:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont --min-read-percent-identity 0.85 -o www/fraction_recovered/long_fraction_recovered.tsv".split()

            run(coverm_cmd)

    if short_reads_2 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -o www/fraction_recovered/short_fraction_recovered.tsv".split()
        run(coverm_cmd)

    elif short_reads_1 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -o www/fraction_recovered/short_fraction_recovered.tsv".split()
        run(coverm_cmd)
        


if __name__ == '__main__':
    os.makedirs('www/fraction_recovered', exist_ok=True)

    long_reads = snakemake.config['long_reads']
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']
    input_fasta = snakemake.input.fasta
    long_read_type = snakemake.config['long_read_type'][0]
    threads = snakemake.threads

    get_fraction_recovered(
        long_reads=long_reads,
        short_reads_1=short_reads_1,
        short_reads_2=short_reads_2,
        input_fasta=input_fasta,
        long_read_type=long_read_type,
        threads=threads,
    )