ruleorder: das_tool_no_vamb > das_tool
ruleorder: das_tool_without_rosella_no_vamb > das_tool_without_rosella
ruleorder: checkm_no_vamb > checkm
ruleorder: checkm_without_rosella_no_vamb > checkm_without_rosella
ruleorder: binner_result_no_vamb > binner_result

# onsuccess:
#     print("Binning finished, no error")
#
# onerror:
#     print("An error occurred")

onstart:
    import os
    import sys


    from snakemake.utils import min_version

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"../../scripts"))

    # minimum required snakemake version
    min_version("6.0")
    long_reads = config["long_reads"]
    fasta = config["fasta"]
    short_reads_1 = config["short_reads_1"]
    short_reads_2 = config["short_reads_2"]
    min_contig_size = config["min_contig_size"]
    min_bin_size = config["min_bin_size"]
    gtdbtk_folder = config["gtdbtk_folder"]
    busco_folder = config["busco_folder"]
    threads = config["max_threads"]
    ## pplacer deadlocks on too many threads
    pplacer_threads = min(48, int(config["pplacer_threads"]))

    if long_reads == "none" and short_reads_1 == "none":
        sys.exit("Need at least one of long_reads or short_reads_1")
    if long_reads != "none" and not os.path.exists(long_reads[0]):
        sys.exit("long_reads does not point to a file")
    if short_reads_1 != "none" and not os.path.exists(short_reads_1[0]):
        sys.exit("short_reads_1 does not point to a file")
    if short_reads_2 != "none" and not os.path.exists(short_reads_2[0]):
        sys.exit("short_reads_2 does not point to a file")
    if gtdbtk_folder != "none" and not os.path.exists(gtdbtk_folder):
        sys.stderr.write("gtdbtk_folder does not point to a folder\n")
    if busco_folder != "none" and not os.path.exists(busco_folder):
        sys.stderr.write("busco_folder does not point to a folder\n")

if config['fasta'] == 'none':
    config['fasta'] = 'assembly/final_contigs.fasta'


rule prepare_binning_files:
    input:
        fasta = config["fasta"]
    output:
        maxbin_coverage = "data/maxbin.cov.list",
        metabat_coverage = "data/coverm.cov",
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    script:
        "../../scripts/get_coverage.py"

rule get_bam_indices:
    input:
        coverage = "data/coverm.cov"
    output:
        bams = "data/binning_bams/done"
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    shell:
        "ls data/binning_bams/*.bam | parallel -j 1 'samtools index -@ {threads} {{}} {{}}.bai' &&"
        "touch data/binning_bams/done"


rule maxbin_binning:
    input:
        fasta = config["fasta"],
        maxbin_cov = "data/maxbin.cov.list"
    params:
        min_contig_size = config["min_contig_size"]
    output:
        "data/maxbin2_bins/done"
    conda:
        "envs/maxbin2.yaml"
    benchmark:
        "benchmarks/maxbin2.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "mkdir -p data/maxbin2_bins && "
        "run_MaxBin.pl -contig {input.fasta} -thread {threads} -abund_list {input.maxbin_cov} "
        "-out data/maxbin2_bins/maxbin -min_contig_length {params.min_contig_size} && "
        "touch data/maxbin2_bins/done"


rule concoct_binning:
    input:
        fasta = config["fasta"],
        bam_done = "data/binning_bams/done"
    params:
        min_contig_size = config["min_contig_size"]
    output:
        "data/concoct_bins/done"
    conda:
        "envs/concoct.yaml"
    benchmark:
        "benchmarks/concoct.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "mkdir -p data/concoct_working && "
        "cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last -b data/concoct_working/contigs_10K.bed > data/concoct_working/contigs_10K.fa && "
        "concoct_coverage_table.py data/concoct_working/contigs_10K.bed data/binning_bams/*.bam > data/concoct_working/coverage_table.tsv && "
        "concoct --threads {threads} -l {params.min_contig_size} --composition_file data/concoct_working/contigs_10K.fa --coverage_file data/concoct_working/coverage_table.tsv -b data/concoct_working/ && "
        "merge_cutup_clustering.py data/concoct_working/clustering_gt{params.min_contig_size}.csv > data/concoct_working/clustering_merged.csv && "
        "mkdir -p data/concoct_bins && "
        "extract_fasta_bins.py {input.fasta} data/concoct_working/clustering_merged.csv --output_path data/concoct_bins/ && "
        "rm -rf data/binning_bams/; "
        "touch data/concoct_bins/done"


rule vamb_jgi_filter:
    input:
        fasta = config["fasta"],
        done = "data/coverm.cov"
    output:
        vamb_bams_done = "data/coverm.filt.cov"
    threads:
        config["max_threads"]
    params:
        min_contig_size = config['min_contig_size']
    run:
        import pandas as pd
        coverm_out = pd.read_csv('data/coverm.cov', sep='\t')
        coverm_out = coverm_out[coverm_out['contigLen']>=int(params.min_contig_size)]
        coverm_out.to_csv("data/coverm.filt.cov", sep='\t', index=False)


rule vamb_binning:
    input:
        coverage = "data/coverm.filt.cov",
        fasta = config["fasta"],
    params:
        min_bin_size = config["min_bin_size"],
        min_contig_size = config["min_contig_size"]
    output:
        "data/vamb_bins/done"
    conda:
        "envs/vamb.yaml"
    benchmark:
        "benchmarks/vamb.benchmark.txt"
    threads:
         config["max_threads"]
    shell:
        "rm -rf data/vamb_bins/; "
        "vamb --outdir data/vamb_bins/ -p {threads} --jgi data/coverm.filt.cov --fasta {input.fasta} --minfasta {params.min_bin_size} -m {params.min_contig_size} || "
        "touch {output[0]}"


rule benchmark_vamb:
    input:
        "data/vamb_bins/done"
    output:
        "data/vamb_bins/checkm.out"
    params:
        pplacer_threads = config["pplacer_threads"],
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna data/vamb_bins/bins/ data/vamb_bins/checkm --tab_table -f data/vamb_bins/checkm.out'


rule vamb_skip:
    output:
        "data/vamb_bins/skipped"
    shell:
        "mkdir -p data/vamb_bins/;"
        "touch data/vamb_bins/skipped"


rule metabat2:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    output:
        metabat_done = "data/metabat_bins_2/done"
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/metabat_2.benchmark.txt"
    shell:
        "metabat -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_2/binned_contigs; "
        "touch data/metabat_bins_2/done"


rule metabat_spec:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    output:
        'data/metabat_bins_spec/done'
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_spec.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --specific -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_spec/binned_contigs; "
        "touch data/metabat_bins_spec/done"

rule metabat_sspec:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    output:
        'data/metabat_bins_sspec/done'
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_sspec.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --superspecific "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sspec/binned_contigs; "
        "touch data/metabat_bins_sspec/done"

rule metabat_sens:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    output:
        'data/metabat_bins_sens/done'
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_sens.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --sensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sens/binned_contigs; "
        "touch data/metabat_bins_sens/done"

rule metabat_ssens:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    output:
        'data/metabat_bins_ssens/done'
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_ssens.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --supersensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_ssens/binned_contigs; "
        "touch data/metabat_bins_ssens/done"

rule rosella:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"]
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    output:
        "data/rosella_bins/done"
    conda:
        "envs/rosella.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/rosella.benchmark.txt"
    shell:
        "rosella bin -r {input.fasta} -i {input.coverage} -t {threads} -o data/rosella_bins "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} --b-tail 0.4; "
        "touch data/rosella_bins/done"



rule das_tool:
    input:
        fasta = config["fasta"],
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        rosella_done = "data/rosella_bins/done",
        vamb_done = "data/vamb_bins/done",
    output:
        das_tool_done = "data/das_tool_bins/done"
    threads:
        config["max_threads"]
    conda:
        "envs/das_tool.yaml"
    benchmark:
        "benchmarks/das_tool.benchmark.txt"
    shell:
        """
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/metabat_bins_2.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/metabat_bins_spec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/vamb_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_bins -e fna > data/rosella_bins.tsv; 
        scaffold2bin_files=$(find data/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g"); 
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} \
         -i $scaffold2bin_files \
         -c {input.fasta} \
         -o data/das_tool_bins/das_tool && \
        touch data/das_tool_bins/done
        """


rule das_tool_no_vamb:
    input:
        fasta = config["fasta"],
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        rosella_done = "data/rosella_bins/done",
        skipped_vamb = "data/vamb_bins/skipped"
    output:
        vamb_skipped = "data/das_tool_bins/skipped_vamb"
    threads:
        config["max_threads"]
    conda:
        "envs/das_tool.yaml"
    shell:
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/metabat_bins_2.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/metabat_bins_spec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/rosella_bins -e fna > data/rosella_bins.tsv && "
        "DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads}"
        " -i data/metabat_bins_2.tsv,data/rosella_bins.tsv,data/metabat_bins_sspec.tsv,data/metabat_bins_spec.tsv,data/metabat_bins_ssens.tsv,data/metabat_bins_sens.tsv,data/maxbin_bins.tsv,data/concoct_bins.tsv"
        " -c {input.fasta} -o data/das_tool_bins/das_tool && "
        "touch data/das_tool_bins/done; "
        "touch data/das_tool_bins/skipped_vamb"


rule get_abundances:
    input:
        "data/singlem_out/singlem_appraise.svg"
    output:
        "data/coverm_abundances.tsv"
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/get_abundances.py"

rule checkm:
    input:
        done = "data/das_tool_bins/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/checkm.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_bins/das_tool_DASTool_bins data/checkm --tab_table -f data/checkm.out'

rule checkm_no_vamb:
    input:
        done = "data/das_tool_bins/skipped_vamb"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/checkm/skipped_vamb"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_bins/das_tool_DASTool_bins data/checkm --tab_table -f data/checkm.out; '
        'touch data/checkm/skipped_vamb; '


rule gtdbtk:
    input:
        done_file = "data/das_tool_bins/done"
    output:
        done = "data/gtdbtk/done"
    params:
        gtdbtk_folder = config['gtdbtk_folder'],
        pplacer_threads = config["pplacer_threads"]        
    conda:
        "../../envs/gtdbtk.yaml"
    threads:
        config["max_threads"]
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_folder} && "
        "gtdbtk classify_wf --cpus {threads} --pplacer_cpus {params.pplacer_threads} --extension fa "
        "--genome_dir data/das_tool_bins/das_tool_DASTool_bins --out_dir data/gtdbtk && touch data/gtdbtk/done"


rule gtdbtk_no_vamb:
    input:
        done_file = "data/das_tool_bins/skipped_vamb"
    output:
        done = "data/gtdbtk/skipped_vamb"
    params:
        gtdbtk_folder = config['gtdbtk_folder'],
        pplacer_threads = config["pplacer_threads"]        
    conda:
        "../../envs/gtdbtk.yaml"
    threads:
        config["max_threads"]
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_folder} && "
        "gtdbtk classify_wf --cpus {threads} --pplacer_cpus {params.pplacer_threads} --extension fa "
        "--genome_dir data/das_tool_bins/das_tool_DASTool_bins --out_dir data/gtdbtk && touch data/gtdbtk/skipped_vamb"


rule binner_result:
    input:
        done = "data/das_tool_without_rosella/done"
    output:
         "data/all_bins/done"
    params:
        pplacer_threads = config['pplacer_threads']
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        "mkdir -p data/all_bins && cd data/all_bins; "
        "ln -s ../vamb_bins/bins/*.fna ./ && ls *.fna | parallel 'mv {{}} vamb_bins_{{}}'; "
        "ln -s ../metabat_bins_2/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat2.fna'; "
        "ln -s ../metabat_bins_sens/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_sens.fna'; "
        "ln -s ../metabat_bins_spec/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_spec.fna'; "
        "ln -s ../metabat_bins_ssens/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_ssens.fna'; "
        "ln -s ../metabat_bins_sspec/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_sspec.fna'; "
        "ln -s ../concoct_bins/*.fa ./ && ls *.fa | parallel 'mv {{}} concoct_{{.}}.fna'; "
        "ln -s ../maxbin2_bins/*.fasta ./ && ls *.fasta | parallel 'mv {{}} maxbin2_{{.}}.fna'; "
        "ln -s ../rosella_bins/*.fna ./; "
        "rm -f \*.fna; "
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna --tab_table ./ checkm > checkm.out; "
        "touch done && cd ../../"

rule binner_result_no_vamb:
    input:
        done = "data/das_tool_without_rosella/skipped_vamb"
    output:
         "data/all_bins/skipped_vamb"
    params:
        pplacer_threads = config['pplacer_threads']
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        "mkdir -p data/all_bins && cd data/all_bins; "
        "ln -s ../metabat_bins_2/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat2.fna'; "
        "ln -s ../metabat_bins_sens/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_sens.fna'; "
        "ln -s ../metabat_bins_spec/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_spec.fna'; "
        "ln -s ../metabat_bins_ssens/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_ssens.fna'; "
        "ln -s ../metabat_bins_sspec/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.metabat_sspec.fna'; "
        "ln -s ../concoct_bins/*.fa ./ && ls *.fa | parallel 'mv {{}} concoct_{{.}}.fna'; "
        "ln -s ../maxbin2_bins/*.fasta ./ && ls *.fasta | parallel 'mv {{}} maxbin2_{{.}}.fna'; "
        "ln -s ../rosella_bins/*.fna ./; "
        "rm -f \*.fna; "
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna --tab_table ./ checkm > checkm.out; "
        "touch skipped_vamb && cd ../../"

rule singlem_pipe_reads:
    output:
        "data/singlem_out/metagenome.combined_otu_table.csv"
    conda:
        "../../envs/singlem.yaml"
    script:
        "../../scripts/singlem_reads.py"

rule singlem_appraise:
    input:
        metagenome = "data/singlem_out/metagenome.combined_otu_table.csv",
        gtdbtk_done = "data/gtdbtk/done"
    output:
        "data/singlem_out/singlem_appraise.svg"
    params:
        pplacer_threads = config['pplacer_threads'],
        fasta = config['fasta']
    conda:
        "../../envs/singlem.yaml"
    shell:
        "singlem pipe --threads {params.pplacer_threads} --sequences data/das_tool_bins/das_tool_DASTool_bins/*.fa --otu_table data/singlem_out/genomes.otu_table.csv; "
        "singlem pipe --threads {params.pplacer_threads} --sequences {params.fasta} --otu_table data/singlem_out/assembly.otu_table.csv; "
        "singlem appraise --metagenome_otu_tables {input.metagenome} --genome_otu_tables data/singlem_out/genomes.otu_table.csv "
        "--assembly_otu_table data/singlem_out/assembly.otu_table.csv "
        "--plot data/singlem_out/singlem_appraise.svg --output_binned_otu_table data/singlem_out/binned.otu_table.csv "
        "--output_unbinned_otu_table data/singlem_out/unbinned.otu_table.csv"

rule singlem_appraise_no_vamb:
    input:
        metagenome = "data/singlem_out/metagenome.combined_otu_table.csv",
        gtdbtk_done = "data/gtdbtk/skipped_vamb"
    output:
        "data/singlem_out/skipped_vamb"
    params:
        pplacer_threads = config['pplacer_threads'],
        fasta = config['fasta']
    conda:
        "../../envs/singlem.yaml"
    shell:
        "singlem pipe --threads {params.pplacer_threads} --sequences data/das_tool_bins/das_tool_DASTool_bins/*.fa --otu_table data/singlem_out/genomes.otu_table.csv; "
        "singlem pipe --threads {params.pplacer_threads} --sequences {params.fasta} --otu_table data/singlem_out/assembly.otu_table.csv; "
        "singlem appraise --metagenome_otu_tables {input.metagenome} --genome_otu_tables data/singlem_out/genomes.otu_table.csv "
        "--assembly_otu_table data/singlem_out/assembly.otu_table.csv "
        "--plot data/singlem_out/singlem_appraise.svg --output_binned_otu_table data/singlem_out/binned.otu_table.csv "
        "--output_unbinned_otu_table data/singlem_out/unbinned.otu_table.csv; "
        "touch data/singlem/skipped_vamb"

rule recover_mags:
    input:
        das = "data/das_tool_bins/done",
        gtdbtk = "data/gtdbtk/done",
        checkm = "data/checkm.out",
        coverm = "data/coverm_abundances.tsv",
        singlem = "data/singlem_out/singlem_appraise.svg"
    conda:
        "../../envs/coverm.yaml"
    output:
        bins = "bins/done",
        taxonomy = "taxonomy/done",
        diversity = 'diversity/done',
        quality = 'bins/checkm.out'
    threads:
        config["max_threads"]
    shell:
        # Use --precluster-method finch so dashing-related install problems are avoided i.e. https://github.com/dnbaker/dashing/issues/41
        "mkdir -p data/pre_galah_bins && cd data/pre_galah_bins/ && rm -f * && ln -s ../das_tool_bins/das_tool_DASTool_bins/* ./ && cd ../../ && "
        "coverm cluster --precluster-method finch -t {threads} --checkm-tab-table data/checkm.out --genome-fasta-directory data/pre_galah_bins/ -x fa --output-representative-fasta-directory data/galah_bins --ani 0.97; "
        "mv data/das_tool_bins/das_tool_DASTool_bins bins/final_bins; "
        "mv data/galah_bins bins/galah_dereplicated; "
        "mv data/checkm.out bins/; "
        "mv data/coverm_abundances.tsv bins/; "
        "mv data/coverm.cov bins/; "
        "mv data/*_bins* bins/; "
        "mv data/singlem_out/ diversity/; "
        "mv data/gtdbtk/ taxonomy/; "
        "touch bins/done; "
        "touch diversity/done; "
        "touch taxonomy/done; "

rule recover_mags_no_vamb:
    input:
        "data/das_tool_bins/skipped_vamb",
        "data/gtdbtk/skipped_vamb",
        "data/checkm/skipped_vamb",
        "data/coverm_abundances.tsv",
        "data/singlem_out/skipped_vamb"
    conda:
        "../../envs/coverm.yaml"
    output:
        "data/done"
    threads:
        config["max_threads"]
    shell:
        # Use --precluster-method finch so dashing-related install problems are avoided i.e. https://github.com/dnbaker/dashing/issues/41
        "mkdir -p data/pre_galah_bins && cd data/pre_galah_bins/ && rm -f * && ln -s ../das_tool_bins/das_tool_DASTool_bins/* ./ && cd ../../ && "
        "coverm cluster --precluster-method finch -t {threads} --checkm-tab-table data/checkm.out "
        "--genome-fasta-directory data/pre_galah_bins/ -x fa --output-representative-fasta-directory data/galah_bins "
        "--ani 0.97 && touch data/done"
