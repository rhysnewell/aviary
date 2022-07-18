
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

rule eggnog:
    input:
        mag_folder = config['mag_directory'],
        # mag_extension = config['mag_extension'],
        eggnog_db = config['eggnog_folder']
    params:
        mag_extension = config['mag_extension']
    group: 'annotation'
    output:
        done = 'data/eggnog/done'
    threads:
        config['max_threads']
    conda:
        'envs/eggnog.yaml'
    shell:
        # 'download_eggnog_data.py --data_dir {input.eggnog_db} -y; '
        'mkdir -p data/eggnog/; '
        'find {input.mag_folder}/*.{params.mag_extension} | parallel -j1 \'emapper.py --data_dir {input.eggnog_db} --dmnd_db {input.eggnog_db}/*dmnd --cpu {threads} -m diamond --itype genome --genepred prodigal -i {{}} --output_dir data/eggnog/ -o {{/.}} || echo "Genome already annotated"\'; '
        'touch data/eggnog/done; '

rule gtdbtk:
    input:
        mag_folder = config['mag_directory']
    group: 'annotation'
    output:
        done = "data/gtdbtk/done"
    params:
        gtdbtk_folder = config['gtdbtk_folder'],
        pplacer_threads = config["pplacer_threads"],
        extension = config['mag_extension']
    conda:
        "../../envs/gtdbtk.yaml"
    threads:
        config["max_threads"]
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_folder} && "
        "gtdbtk classify_wf --cpus {threads} --pplacer_cpus {params.pplacer_threads} --extension {params.extension} "
        "--genome_dir {input.mag_folder} --out_dir data/gtdbtk && touch data/gtdbtk/done"

rule annotate:
    input:
         'data/gtdbtk/done',
         'data/eggnog/done',
    output:
         'annotation/done',
    shell:
         """
         ln -sr data/gtdbtk taxonomy; 
         ln -sr data/eggnog annotation; 
         touch annotation/done;
         """
