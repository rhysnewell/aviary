# onsuccess:
#     print("Strain analysis finished, no error")
#
# onerror:
#     print("An error occurred")


onstart:
    import os
    import sys

    from snakemake.utils import logger, min_version

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "../../scripts"))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))

    # minimum required snakemake version
    min_version("6.0")
    long_reads = config["long_reads"]
    fasta = config["fasta"]
    short_reads_1 = config["short_reads_1"]
    short_reads_2 = config["short_reads_2"]
    min_contig_size = config["min_contig_size"]
    min_bin_size = config["min_bin_size"]
    gtdbtk_folder = config["gtdbtk_folder"]
    busco_folder = config['busco_folder']
    threads = config["max_threads"]
    ## pplacer deadlocks on too many threads
    pplacer_threads = min(48, int(config["pplacer_threads"]))

    if gtdbtk_folder != "none" and not os.path.exists(gtdbtk_folder):
        sys.stderr.write("gtdbtk_folder does not point to a folder\n")
    if busco_folder != "none" and not os.path.exists(busco_folder):
        sys.stderr.write("busco_folder does not point to a folder\n")


if config['fasta'] == 'none':
    config['fasta'] = 'assembly/final_contigs.fasta'

if config['mag_directory'] == 'none':
    config['mag_directory'] = 'bins/final_bins'
if config['mag_extension'] == 'none':
    config['mag_extension'] = 'fna'

# rule generate_bams:
#     input:
#         bins_directory = 'data/galah_bins/',
#     output:
#         finished_mapping = 'data/binned_bams/{sample}.bam'
#     threads:
#         config['max_threads']
#     conda:
#         '../../envs/coverm.yaml'
#     shell:
#         'coverm genome '


rule lorikeet:
    input:
         mag_directory = config['mag_directory'],
    output:
         output_directory = directory("strain_diversity/")
    conda:
        "envs/lorikeet.yaml"
    params:
        mag_extension = config['mag_extension'],
        parallel_genomes = 8,
    resources:
        mem_mb=int(config["max_memory"])*512
    threads:
        config["max_threads"]
    script:
        "scripts/run_lorikeet.py"