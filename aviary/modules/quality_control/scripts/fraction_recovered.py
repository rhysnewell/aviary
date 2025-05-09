from subprocess import run, STDOUT
import os
import argparse


def get_fraction_recovered(
    long_reads,
    short_reads_1,
    short_reads_2,
    input_fasta: str,
    long_read_type: str,
    threads: int,
    log: str,
):
    if long_reads != "none" and not os.path.exists("data/long_cov.tsv"):
        if long_read_type in ["ont", "ont_hq"]:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont --min-read-percent-identity 0.85 -o www/fraction_recovered/long_fraction_recovered.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

        elif long_read_type in ["rs", "sq", "ccs", "hifi"]:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-pb --min-read-percent-identity 0.9 -o www/fraction_recovered/long_fraction_recovered.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

        else:
            coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont --min-read-percent-identity 0.85 -o www/fraction_recovered/long_fraction_recovered.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT)

    if short_reads_2 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -o www/fraction_recovered/short_fraction_recovered.tsv".split()

        with open(log, "a") as logf:
            run(coverm_cmd, stdout=logf, stderr=STDOUT)

    elif short_reads_1 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = f"coverm genome -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -o www/fraction_recovered/short_fraction_recovered.tsv".split()

        with open(log, "a") as logf:
            run(coverm_cmd, stdout=logf, stderr=STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate fraction recovered from sequencing reads')
    parser.add_argument('--input-fasta', required=True, help='Input FASTA file')
    parser.add_argument('--long-reads', nargs='+', default=['none'], help='Long read files')
    parser.add_argument('--short-reads-1', nargs='+', default=['none'], help='Short reads (first pair or interleaved)')
    parser.add_argument('--short-reads-2', nargs='+', default=['none'], help='Short reads (second pair)')
    parser.add_argument('--long-read-type', default='ont', help='Long read type (ont, ont_hq, rs, sq, ccs, hifi)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--log', default='get_fraction_recovered.log', help='Log file')
    
    args = parser.parse_args()
    
    os.makedirs('www/fraction_recovered', exist_ok=True)

    with open(args.log, "w") as logf: 
        pass

    get_fraction_recovered(
        long_reads=args.long_reads,
        short_reads_1=args.short_reads_1,
        short_reads_2=args.short_reads_2,
        input_fasta=args.input_fasta,
        long_read_type=args.long_read_type,
        threads=args.threads,
        log=args.log,
    )