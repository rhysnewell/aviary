from subprocess import run

def run_flye(
    long_read_type: str,
    input_fastq: str,
    output_dir: str,
    meta_flag: bool,
    threads: int
):
    meta = ""
    flye_type = "--nano-raw"
    if meta_flag:
        meta = "--meta"
    if long_read_type == 'ont':
        flye_type = "--nano-raw"
    elif long_read_type == 'ont_hq':
        flye_type = "--nano-hq"
    elif long_read_type == 'ccs' or long_read_type == 'hifi':
        flye_type = "--pacbio-hifi"
    else:
        flye_type = "--pacbio-raw"
    
    flye_cmd = f"flye {flye_type} {input_fastq} {meta} -o {output_dir} -t {threads}".split()
    run(flye_cmd)


if __name__ == '__main__':
    long_read_type = snakemake.params.long_read_type
    input_fastq = snakemake.input.fastq
    output_dir = "data/flye"
    meta_flag = True
    threads = snakemake.threads

    run_flye(
        long_read_type=long_read_type,
        input_fastq=input_fastq,
        output_dir=output_dir,
        meta_flag=meta_flag,
        threads=threads
    )