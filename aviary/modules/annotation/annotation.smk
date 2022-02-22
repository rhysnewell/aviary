
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

rule prodigal:
    input:
        fasta = config['fasta']
    output:
        "data/genes.gff"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal -i {input.fasta} -f gff -o {output} -p meta"


rule busco_mag_directory:
    input:
        mags = config['mag_directory']
    output:
        done = 'data/busco/done'
    params:
        busco_folder = config["busco_folder"],
        ext = config['mag_extension']
    conda:
        "envs/busco.yaml"
    threads:
        config["max_threads"]
    script:
        """
        mkdir -p data/busco && \
        cd data/busco && \
        minimumsize=500000 && \
        for file in {input.mags}/*{params.ext};do 
            actualsize=$(wc -c <\"$file\"); 
            if [ $actualsize -ge $minimumsize ]; then 
                if [ ! -d bacteria_odb10.${{file:39:-3}} ]; then
                    busco -q -c {threads} -i $file -o bacteria_odb10.${{file:39:-3}} \
                    -l {params.busco_folder}/bacteria_odb10 -m geno;
                fi
                if [ ! -d eukaryota_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o eukaryota_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/eukaryota_odb10 -m geno; 
                fi
                if [ ! -d embryophyta_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o embryophyta_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/embryophyta_odb10 -m geno; 
                fi
                if [ ! -d fungi_odb10.${{file:39:-3}} ]; then
                busco -q -c {threads} -i $file -o fungi_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/fungi_odb10 -m geno; 
                fi 
                # busco -q -c {threads} -i $file -o metazoa_odb10.${{file:39:-3}} \
                -l {params.busco_folder}/metazoa_odb10 -m geno; 
                # busco -q -c {threads} -i $file -o protists_ensembl.${{file:39:-3}} \
                -l {params.busco_folder}/protists_ensembl -m geno; 
            fi
        done && \
        cd ../../ && \
        touch data/busco/done
        """

# rule enrichm_directory:
#     input:
#         mags = config['mag_directory']
#     output:
#         out_folder = "data/enrichm",
#         done = "data/enrichm/done"
#     params:
#         busco_folder = config["enrichm_folder"],
#         ext = config['mag_extension']
#     conda:
#         "envs/enrichm.yaml"
#     threads:
#         config["max_threads"]
#     shell:
#         """
#         enrichm annotate --output {output.out_folder} --genome_directory {input.mags} --suffix {params.ext} --ko_hmm --pfam --tigrfam --clusters --orthogroup --cazy --ec --parallel {threads} &&
#         enrichm classify --output {output.out_folder} --genome_and_annotation_matrix {output.out_folder}/ko_frequency_table.tsv;
#         touch data/enrichm/done
#         """

rule gtdbtk:
    input:
        mag_folder = config['mag_directory']
    group: 'binning'
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

rule complete_annotation:
    input:
         'data/gtdbtk/done',
         'data/enrichm/done',
    output:
         'annotation/done',
    shell:
         """
         mkdir -p annotation;
         ln -s data/enrichm annotation/enrichm
         ln -s data/busco annotation/busco
         """