ruleorder: busco_bins_provided > busco

# onsuccess:
#     print("Annotation finished, no error")
#
# onerror:
#     print("An error occurred")

onstart:
    import os
    import sys

    from snakemake.utils import logger, min_version

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "../../scripts"))

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

rule prodigal:
    input:
        fasta = config['fasta']
    output:
        "data/genes.gff"
    conda:
        "../../envs/prodigal.yaml"
    shell:
        "prodigal -i {input.fasta} -f gff -o {output} -p meta"

# Run BUSCO on the bins
rule busco_bins_provided:
    input:
        mags = config['mags'],
        # ext = config
    output:
        done = "data/busco/done"
    params:
        busco_folder = config["busco_folder"]
    conda:
        "../../envs/busco.yaml"
    threads:
        config["max_threads"]
    script:
        """
        mkdir -p data/busco && \
        cd data/busco && \
        minimumsize=500000 && \
        for file in ../das_tool_bins/das_tool_DASTool_bins/*.fa;do 
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

# Run BUSCO on the bins
rule busco:
    input:
        mags = "bins/done",

    output:
        done = "data/busco/done"
    params:
        busco_folder = config["busco_folder"]
    conda:
        "../../envs/busco.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        mkdir -p data/busco && \
        cd data/busco && \
        minimumsize=500000 && \
        for file in ../../bins/final_bins/*.fa;do 
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
