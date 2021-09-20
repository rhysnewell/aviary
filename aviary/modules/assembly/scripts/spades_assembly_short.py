import subprocess

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen(
        "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -1 %s -2 %s" %
        (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0],
         snakemake.config["short_reads_2"][0]), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    subprocess.Popen(
        "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -1 %s" %
        (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0]),
        shell=True).wait()

# subprocess.Popen(
#     "ln -s data/spades_assembly/scaffolds.fasta data/final_contigs.fasta", shell=True
# ).wait()
