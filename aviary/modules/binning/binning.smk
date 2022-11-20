ruleorder: dereplicate_and_get_abundances_paired > dereplicate_and_get_abundances_interleaved
ruleorder: checkm_rosella > amber_checkm_output
ruleorder: checkm_metabat2 > amber_checkm_output
ruleorder: checkm_semibin > amber_checkm_output

onstart:
    from snakemake.utils import min_version

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"../../scripts"))
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

import os
import sys
import glob

rule prepare_binning_files:
    input:
        fasta = config["fasta"]
    group: 'binning'
    output:
        maxbin_coverage = "data/maxbin.cov.list",
        metabat_coverage = "data/coverm.cov"
    params:
        tmpdir = config['tmpdir']
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/get_coverage.py"


rule get_bam_indices:
    input:
        coverage = "data/coverm.cov"
    group: 'binning'
    output:
        bams = "data/binning_bams/done"
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    shell:
        "ls data/binning_bams/*.bam | parallel -j 1 'samtools index -@ {threads} {{}} {{}}.bai' &&"
        "touch data/binning_bams/done"


rule maxbin2:
    input:
        fasta = ancient(config["fasta"]),
        maxbin_cov = ancient("data/maxbin.cov.list")
    params:
        min_contig_size = config["min_contig_size"]
    resources:
        mem_mb=int(config["max_memory"])*128
    group: 'binning'
    output:
        "data/maxbin2_bins/done"
    conda:
        "envs/maxbin2.yaml"
    benchmark:
        "benchmarks/maxbin2.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/maxbin2_bins/; "
        "mkdir -p data/maxbin2_bins && "
        "run_MaxBin.pl -contig {input.fasta} -thread {threads} -abund_list {input.maxbin_cov} "
        "-out data/maxbin2_bins/maxbin -min_contig_length {params.min_contig_size} && "
        "touch {output[0]} || touch {output[0]}"


rule concoct:
    input:
        fasta = ancient(config["fasta"]),
        bam_done = ancient("data/binning_bams/done")
    params:
        min_contig_size = config["min_contig_size"]
    resources:
        mem_mb=int(config["max_memory"])*128
    group: 'binning'
    output:
        "data/concoct_bins/done"
    conda:
        "envs/concoct.yaml"
    benchmark:
        "benchmarks/concoct.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/concoct_*/; "
        "mkdir -p data/concoct_working && "
        "cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last -b data/concoct_working/contigs_10K.bed > data/concoct_working/contigs_10K.fa && "
        "concoct_coverage_table.py data/concoct_working/contigs_10K.bed data/binning_bams/*.bam > data/concoct_working/coverage_table.tsv && "
        "concoct --threads {threads} -l {params.min_contig_size} --composition_file data/concoct_working/contigs_10K.fa --coverage_file data/concoct_working/coverage_table.tsv -b data/concoct_working/ 2>/dev/null && "
        "merge_cutup_clustering.py data/concoct_working/clustering_gt{params.min_contig_size}.csv > data/concoct_working/clustering_merged.csv && "
        "mkdir -p data/concoct_bins && "
        "extract_fasta_bins.py {input.fasta} data/concoct_working/clustering_merged.csv --output_path data/concoct_bins/ && "
        "touch {output[0]} || touch {output[0]}"


rule vamb_jgi_filter:
    """
    vamb has to have to coverage file filtered prior to running otherwise it throws an error
    Outputs a coverage file containing no contigs smaller than minimum contig size
    """
    input:
        fasta = ancient(config["fasta"]),
        done = ancient("data/coverm.cov")
    group: 'binning'
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


rule vamb:
    """
    Perform binning via vamb. Vamb frequently breaks and won't produce any bins or errors out. As such, whenenver 
    vamb throws an error this rule will catch it and create the output regardless. You'll know if vamb failed as there
    will be no bins produced by it but all other files will be there
    """
    input:
        coverage = ancient("data/coverm.filt.cov"),
        fasta = ancient(config["fasta"]),
    params:
        min_bin_size = config["min_bin_size"],
        min_contig_size = config["min_contig_size"],
        vamb_threads = int(config["max_threads"]) // 2 # vamb use double the threads you give it
    resources:
        mem_mb=int(config["max_memory"])*128
    group: 'binning'
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
        "bash -c 'vamb --outdir data/vamb_bins/ -p {params.vamb_threads} --jgi {input.coverage} --fasta {input.fasta} "
        "--minfasta {params.min_bin_size} -m {params.min_contig_size} && touch {output[0]}' || "
        "touch {output[0]} && mkdir -p data/vamb_bins/bins"


rule vamb_skip:
    group: 'binning'
    output:
        "data/vamb_bins/skipped"
    shell:
        "mkdir -p data/vamb_bins/;"
        "touch data/vamb_bins/skipped"


rule metabat2:
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"]
    group: 'binning'
    output:
        metabat_done = "data/metabat_bins_2/done"
    conda:
        "envs/metabat2.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/metabat_2.benchmark.txt"
    shell:
        "rm -rf data/metabat_bins_2/; "
        "metabat -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_2/binned_contigs && "
        "touch {output[0]} || touch {output[0]}"


rule metabat_spec:
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    group: 'binning'
    output:
        'data/metabat_bins_spec/done'
    conda:
        "envs/metabat2.yaml"
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_spec.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/metabat_bins_spec; "
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --specific -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_spec/binned_contigs && "
        "touch {output[0]} || touch {output[0]}"

rule metabat_sspec:
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    group: 'binning'
    output:
        'data/metabat_bins_sspec/done'
    conda:
        "envs/metabat2.yaml"
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_sspec.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/metabat_bins_sspec; "
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --superspecific "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sspec/binned_contigs && "
        "touch {output[0]} || touch {output[0]}"

rule metabat_sens:
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    group: 'binning'
    output:
        'data/metabat_bins_sens/done'
    conda:
        "envs/metabat2.yaml"
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_sens.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/metabat_bins_sens; "
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --sensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sens/binned_contigs && "
        "touch {output[0]} || touch {output[0]}"

rule metabat_ssens:
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    group: 'binning'
    output:
        'data/metabat_bins_ssens/done'
    conda:
        "envs/metabat2.yaml"
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"]
    benchmark:
        "benchmarks/metabat_ssens.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "rm -rf data/metabat_bins_ssens; "
        "metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --supersensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_ssens/binned_contigs && "
        "touch {output[0]} || touch {output[0]}"

rule rosella:
    """
    Runs Rosella.
    """
    input:
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"])
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    resources:
        mem_mb=int(config["max_memory"])*128
    group: 'binning'
    output:
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv",
        done = "data/rosella_bins/done"
    conda:
        "envs/rosella.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/rosella.benchmark.txt"
    shell:
        "rm -rf data/rosella_bins/; "
        "rosella recover -r {input.fasta} -i {input.coverage} -t {threads} -o data/rosella_bins "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} --n-neighbors 200 && "
        "touch {output.done} || touch {output.done}"


rule semibin:
    input:
        fasta = ancient(config["fasta"]),
        bams_indexed = ancient("data/binning_bams/done")
    group: 'binning'
    params:
        # Can't use premade model with multiple samples, so disregard if provided
        semibin_model = config['semibin_model']
    resources:
        mem_mb=int(config["max_memory"])*1024
    output:
        done = "data/semibin_bins/done"
    threads:
        config["max_threads"]
    conda:
        "envs/semibin.yaml"
    benchmark:
        "benchmarks/semibin.benchmark.txt"
    shell:
        "rm -rf data/semibin_bins/; "
        "mkdir -p data/semibin_bins/output_recluster_bins/; "
        "SemiBin single_easy_bin -i {input.fasta} -b data/binning_bams/*.bam -o data/semibin_bins --environment {params.semibin_model} -p {threads} && "
        "touch {output.done} || SemiBin single_easy_bin -i {input.fasta} -b data/binning_bams/*.bam -o data/semibin_bins -p {threads} "
        "&& touch {output.done} || touch {output.done}"

rule checkm_rosella:
    input:
        done = ancient("data/rosella_bins/done")
    params:
        pplacer_threads = config["pplacer_threads"],
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/rosella_bins/",
        extension = "fna"
    group: 'binning'
    output:
        output_folder = directory("data/rosella_bins/checkm2_out/"),
        output_file = "data/rosella_bins/checkm.out"
    conda:
        "../../envs/checkm2.yaml"
    threads:
        config["max_threads"]
    shell:
        'touch {output.output_file}; '
        'if [ `ls "{params.bin_folder}" |grep .fna$ |wc -l` -eq 0 ]; then '
        'echo "No bins found in {params.bin_folder}"; '
        'touch {output.output_file}; '
        'mkdir -p {output.output_folder}; '
        'else '

        'export CHECKM2DB={params.checkm2_db_path}/uniref100.KO.1.dmnd; '
        'echo "Using CheckM2 database $CHECKM2DB"; '
        'checkm2 predict -i {params.bin_folder}/ -x {params.extension} -o {output.output_folder} -t {threads} --force; '
        'cp {output.output_folder}/quality_report.tsv {output.output_file}; '

        'fi'

rule checkm_metabat2:
    input:
        done = ancient("data/metabat_bins_2/done")
    params:
        pplacer_threads = config["pplacer_threads"],
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/metabat_bins_2/",
        extension = "fa"
    group: 'binning'
    output:
        output_folder = directory("data/metabat_bins_2/checkm2_out/"),
        output_file = "data/metabat_bins_2/checkm.out"
    conda:
        "../../envs/checkm2.yaml"
    threads:
        config["max_threads"]
    shell:
        'touch {output.output_file}; '
        'export CHECKM2DB={params.checkm2_db_path}/uniref100.KO.1.dmnd; '
        'echo "Using CheckM2 database $CHECKM2DB"; '
        'checkm2 predict -i {params.bin_folder}/ -x {params.extension} -o {output.output_folder} -t {threads} --force; '
        'cp {output.output_folder}/quality_report.tsv {output.output_file}'

rule checkm_semibin:
    input:
        done = ancient("data/semibin_bins/done")
    params:
        pplacer_threads = config["pplacer_threads"],
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/semibin_bins/output_recluster_bins/",
        extension = "fa"
    group: 'binning'
    output:
        output_folder = directory("data/semibin_bins/checkm2_out/"),
        output_file = "data/semibin_bins/checkm.out"
    conda:
        "../../envs/checkm2.yaml"
    threads:
        config["max_threads"]
    shell:
        'touch {output.output_file}; '
        'export CHECKM2DB={params.checkm2_db_path}/uniref100.KO.1.dmnd; '
        'echo "Using CheckM2 database $CHECKM2DB"; '
        'checkm2 predict -i {params.bin_folder}/ -x {params.extension} -o {output.output_folder} -t {threads} --force; '
        'cp {output.output_folder}/quality_report.tsv {output.output_file}'

rule refine_rosella:
    input:
        checkm = ancient('data/rosella_bins/checkm.out'),
        rosella = ancient('data/rosella_bins/done'),
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"]),
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/rosella_refined/done'
    benchmark:
        'benchmarks/refine_rosella.benchmark.txt'
    params:
        bin_folder = "data/rosella_bins/",
        extension = "fna",
        output_folder = "data/rosella_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 5,
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 15,
        final_refining = False
    resources:
        mem_mb=int(config["max_memory"])
    threads:
        config["max_threads"]
    conda:
        "envs/rosella.yaml"
    script:
        "scripts/rosella_refine.py"

rule refine_metabat2:
    input:
        checkm = ancient('data/metabat_bins_2/checkm.out'),
        rosella = ancient('data/metabat_bins_2/done'),
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"]),
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/metabat2_refined/done'
    resources:
        mem_mb=int(config["max_memory"])
    benchmark:
        'benchmarks/refine_metabat2.benchmark.txt'
    params:
        bin_folder = "data/metabat_bins_2/",
        extension = "fa",
        output_folder = "data/metabat2_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 5,
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 15,
        final_refining = False
    threads:
        config["max_threads"]
    conda:
        "envs/rosella.yaml"
    script:
        "scripts/rosella_refine.py"

rule refine_semibin:
    input:
        checkm = ancient('data/semibin_bins/checkm.out'),
        rosella = ancient('data/semibin_bins/done'),
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"]),
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    resources:
        mem_mb=int(config["max_memory"])
    output:
        'data/semibin_refined/done'
    benchmark:
        'benchmarks/refine_semibin.benchmark.txt'
    params:
        bin_folder = "data/semibin_bins/output_recluster_bins/",
        extension = "fa",
        output_folder = "data/semibin_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 5,
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 15,
        final_refining = False
    threads:
        config["max_threads"]
    conda:
        "envs/rosella.yaml"
    script:
        "scripts/rosella_refine.py"

rule amber_checkm_output:
    input:
        amber_done = "data/amber_refine/for_refine/index.html"
    output:
        metabat_checkm = 'data/metabat_bins_2/checkm.out',
        rosella_checkm = 'data/rosella_bins/checkm.out',
        semibin_checkm = 'data/semibin_bins/checkm.out'
        # dastool_checkm = 'data/das_tool_bins_with_refine/checkm.out'
    run:
        import pandas as pd

        def amber_to_checkm_like(amber_table, prefix="data/rosella_bins/", output="data/rosella_bins/checkm.out", extension="fna"):
            amber_table["Bin Id"] = amber_table["Bin ID"].replace(prefix, "", regex=True)
            amber_table["Bin Id"] = amber_table["Bin Id"].replace(f".{extension}", "", regex=True)
            amber_table["Completeness"] = amber_table["Completeness (bp)"] * 100
            amber_table["Contamination"] = round((1 - amber_table["Purity (bp)"]) * 100, 2)
            amber_table.to_csv(output, sep='\t', index=False)

        try:
            rosella_amber = pd.read_csv("data/amber_refine/for_refine/genome/rosella_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(rosella_amber, "data/rosella_bins/", "data/rosella_bins/checkm.out", "fna")
        except FileNotFoundError:
            pass

        try:
            metabat_amber = pd.read_csv("data/amber_refine/for_refine/genome/m2_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(metabat_amber, "data/metabat_bins_2/", "data/metabat_bins_2/checkm.out", "fa")
        except FileNotFoundError:
            pass

        try:
            semibin_amber = pd.read_csv("data/amber_refine/for_refine/genome/semibin_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(semibin_amber, "data/semibin_bins/output_recluster_bins/", "data/semibin_bins/checkm.out", "fa")
        except FileNotFoundError:
            pass


rule das_tool:
    """
    Runs dasTool on the output of all binning algorithms. If a binner failed to produce bins then their output is ignored
    """
    input:
        fasta = ancient(config["fasta"]),
        metabat2_done = "data/metabat2_refined/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        rosella_done = "data/rosella_refined/done",
        semibin_done = "data/semibin_refined/done",
        vamb_done = "data/vamb_bins/done",
    resources:
        mem_mb=int(config["max_memory"])
    group: 'binning'
    output:
        das_tool_done = "data/das_tool_bins_pre_refine/done"
    threads:
        config["max_threads"]
    conda:
        "envs/das_tool.yaml"
    benchmark:
        "benchmarks/das_tool.benchmark.txt"
    shell:
        """
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/metabat_bins_spec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/vamb_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_refined/final_bins/ -e fna > data/rosella_refined_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat2_refined/final_bins/ -e fna > data/metabat2_refined_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/semibin_refined/final_bins/ -e fna > data/semibin_refined_bins.tsv; 
        scaffold2bin_files=$(find data/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g"); 
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} --score_threshold -42 \
         -i $scaffold2bin_files \
         -c {input.fasta} \
         -o data/das_tool_bins_pre_refine/das_tool && \
        touch data/das_tool_bins_pre_refine/done
        """

rule refine_dastool:
    input:
        checkm = 'data/das_tool_bins_pre_refine/checkm.out',
        das_tool = 'data/das_tool_bins_pre_refine/done',
        coverage = ancient("data/coverm.cov"),
        fasta = ancient(config["fasta"]),
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    resources:
        mem_mb=int(config["max_memory"])
    output:
        temporary('bins/checkm.out'),
        directory('bins/final_bins')
    benchmark:
        'benchmarks/refine_dastool.benchmark.txt'
    params:
        bin_folder = "data/das_tool_bins_pre_refine/das_tool_DASTool_bins/",
        extension = "fa",
        output_folder = "data/refined_bins/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 5,
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 15,
        final_refining = True
    threads:
        config["max_threads"]
    conda:
        "envs/rosella.yaml"
    script:
        "scripts/rosella_refine.py"

rule get_abundances:
    input:
        "bins/checkm.out"
    group: 'binning'
    resources:
        mem_mb=int(config["max_memory"])*1024
    output:
        "data/coverm_abundances.tsv"
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/get_abundances.py"

rule finalize_stats:
    input:
        checkm1_done = "bins/checkm.out",
        checkm2_done = "bins/checkm2_output/quality_report.tsv",
        coverage_file = "data/coverm_abundances.tsv",
        gtdbtk_done = "data/gtdbtk/done"
    output:
        bin_stats = "bins/bin_info.tsv",
        checkm_minimal = "bins/checkm_minimal.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO
        import os

        def find_circular(checkm_output, checkm1=True):
            if checkm1:
                bin_column = "Bin Id"
            else:
                bin_column = "Name"

            circular_contigs = []
            circular_bps = []
            circular_fractions = []

            assembly_info = pd.read_csv("data/flye/assembly_info.txt", sep="\t")

            for bin_name in checkm_output[bin_column]:

                fasta_path = f"bins/final_bins/{bin_name}.fna"
                circular = 0
                circular_bases = 0
                total_size = 0
                for sequence in SeqIO.parse(open(fasta_path), "fasta"):
                    total_size += len(sequence.seq)
                    seq_name = sequence.id.strip("_pilon")

                    if seq_name not in assembly_info["#seq_name"].values:
                        continue

                    found = assembly_info[assembly_info["#seq_name"] == seq_name]

                    if found["circ."].values[0] == "N":
                        continue
                    circular += 1
                    circular_bases += found["length"].values[0]

                circular_bps.append(circular_bases)
                circular_contigs.append(circular)
                circular_fractions.append(circular_bases / total_size)

            checkm_output["Circular contigs"], checkm_output["Circular bp"], checkm_output["Circular fraction"] = [circular_contigs, circular_bps, circular_fractions]
            return checkm_output

        def get_taxonomy(rename_columns="Bin Id"):
            taxa = []
            try:
                df_bac = pd.read_csv(glob.glob("data/gtdbtk/gtdbtk.bac*.summary.tsv")[0], sep="\t")
                taxa.append(df_bac)
            except (FileNotFoundError, IndexError):
                pass

            try:
                df_arc = pd.read_csv(glob.glob("data/gtdbtk/gtdbtk.ar*.summary.tsv")[0], sep="\t")
                taxa.append(df_arc)
            except (FileNotFoundError, IndexError) as e:
                pass

            taxa = pd.concat(taxa)
            taxa.rename({'user_genome' : rename_columns}, inplace=True, axis=1)
            return taxa

        coverage_file = pd.read_csv(input.coverage_file, sep='\t')

        # checkm file for all bins
        checkm1_output = pd.read_csv(input.checkm1_done, sep='\t', comment="[")

        checkm2_output = pd.read_csv(input.checkm2_done, sep='\t')

        checkm1_output.rename({'Completeness' : 'Completeness (CheckM1)', 'Contamination' : 'Contamination (CheckM1)'}, inplace=True, axis=1)
        checkm2_output.rename({'Name' : checkm1_output.columns[0], 'Completeness' : 'Completeness (CheckM2)', 'Contamination' : 'Contamination (CheckM2)'}, inplace=True, axis=1)

        checkm_output = pd.merge(checkm1_output, checkm2_output, on=[checkm1_output.columns[0]])
        is_checkm1 = "Bin Id" in checkm_output.columns
        coverage_file.rename({"Genome" : checkm_output.columns[0]}, inplace=True, axis=1)


        if os.path.isfile("data/flye/assembly_info.txt"):
            checkm_output = find_circular(checkm_output, is_checkm1)

        taxa = get_taxonomy(checkm_output.columns[0])

        merged_out = pd.merge(checkm_output, coverage_file, on=[checkm_output.columns[0]])
        merged_out = pd.merge(merged_out, taxa, on=[checkm_output.columns[0]])
        merged_out.to_csv(output.bin_stats, sep='\t', index=False)

        checkm_minimal = checkm_output[["Bin Id",  "Marker lineage",  "# genomes", "# markers", "# marker sets",
                                        "0", "1", "2", "3", "4", "5+", "Completeness (CheckM1)", "Contamination (CheckM1)",
                                        "Completeness (CheckM2)", "Contamination (CheckM2)", "Strain heterogeneity"]]

        checkm_minimal.to_csv(output.checkm_minimal, sep="\t", index=False)



rule checkm_das_tool:
    input:
        done = "data/das_tool_bins_pre_refine/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    group: 'binning'
    output:
        "data/das_tool_bins_pre_refine/checkm.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_bins_pre_refine/das_tool_DASTool_bins data/das_tool_bins_pre_refine/checkm --tab_table '
        '-f data/das_tool_bins_pre_refine/checkm.out; '
        'checkm qa -o 2 --tab_table -f data/das_tool_bins_pre_refine/checkm.out '
        'data/das_tool_bins_pre_refine/checkm/lineage.ms data/das_tool_bins_pre_refine/checkm/; '


rule singlem_pipe_reads:
    group: 'binning'
    output:
        "data/singlem_out/metagenome.combined_otu_table.csv"
    conda:
        "../../envs/singlem.yaml"
    script:
        "../../scripts/singlem_reads.py"

rule singlem_appraise:
    input:
        metagenome = "data/singlem_out/metagenome.combined_otu_table.csv",
        gtdbtk_done = "data/gtdbtk/done",
        bins_complete = "bins/checkm.out"
    group: 'binning'
    output:
        "data/singlem_out/singlem_appraisal.tsv"
    params:
        pplacer_threads = config['pplacer_threads'],
        fasta = config['fasta']
    threads:
        config["pplacer_threads"]
    conda:
        "../../envs/singlem.yaml"
    log:
        "data/singlem_out/singlem_log.txt"
    shell:
        "singlem pipe --threads {threads} --sequences bins/final_bins/*.fna --otu_table data/singlem_out/genomes.otu_table.csv; "
        "singlem pipe --threads {threads} --sequences {params.fasta} --otu_table data/singlem_out/assembly.otu_table.csv; "
        "singlem appraise --metagenome_otu_tables {input.metagenome} --genome_otu_tables data/singlem_out/genomes.otu_table.csv "
        "--assembly_otu_table data/singlem_out/assembly.otu_table.csv "
        "--plot data/singlem_out/singlem_appraise.svg --output_binned_otu_table data/singlem_out/binned.otu_table.csv "
        "--output_unbinned_otu_table data/singlem_out/unbinned.otu_table.csv 1> data/singlem_out/singlem_appraisal.tsv 2> {log} || "
        "echo 'SingleM Errored, please check data/singlem_out/singlem_log.txt' && touch data/singlem_out/singlem_appraisal.tsv"


rule recover_mags:
    input:
        final_bins = "bins/bin_info.tsv",
        gtdbtk = "data/gtdbtk/done",
        coverm = "data/coverm_abundances.tsv",
        singlem = "data/singlem_out/singlem_appraisal.tsv"
    conda:
        "../../envs/coverm.yaml"
    group: 'binning'
    output:
        bins = "bins/done",
        diversity = 'diversity/done'
    threads:
        config["max_threads"]
    shell:
        "cd bins/; "
        "ln -s ../data/coverm_abundances.tsv ./; "
        "ln -s ../data/coverm.cov ./; "
        "cd ../; "
        "ln -sr data/singlem_out/ diversity || echo 'SingleM linked'; "
        "ln -sr data/gtdbtk taxonomy || echo 'GTDB-tk linked'; "
        "touch bins/done; "
        "touch diversity/done; "
        "rm -f data/binning_bams/*bam; "
        "rm -f data/binning_bams/*bai; "

rule recover_mags_no_singlem:
    input:
        final_bins = "bins/bin_info.tsv",
        coverm = "data/coverm_abundances.tsv",
    conda:
        "../../envs/coverm.yaml"
    group: 'binning'
    output:
        bins = "bins/done",
    threads:
        config["max_threads"]
    shell:
        "cd bins/; "
        "ln -s ../data/coverm_abundances.tsv ./; "
        "ln -s ../data/coverm.cov ./; "
        "cd ../; "
        "touch bins/done; "
        "rm -f data/binning_bams/*bam; "
        "rm -f data/binning_bams/*bai; "

# Special rule to help out with a buggy output
rule dereplicate_and_get_abundances_paired:
    input:
        pe_1 = config["short_reads_1"],
        pe_2 = config["short_reads_2"]
    output:
        output_abundances = 'bins/coverm_abundances.tsv'
    params:
        final_bins = 'bins/final_bins',
        derep_ani = 0.97
    threads:
        config['max_threads']
    conda:
        "../../envs/coverm.yaml"
    shell:
        "coverm genome -t {threads} -d bins/final_bins/ -1 {input.pe_1} -2 {input.pe_2} --min-covered-fraction 0.0 -x fna > bins/coverm_abundances.tsv; "

# Special rule to help out with a buggy output
rule dereplicate_and_get_abundances_interleaved:
    input:
        pe_1 = config["short_reads_1"],
    output:
        output_abundances = 'bins/coverm_abundances.tsv'
    params:
        final_bins = 'bins/final_bins',
        derep_ani = 0.97
    threads:
        config['max_threads']
    conda:
        "../../envs/coverm.yaml"
    shell:
        "coverm genome -t {threads} -d bins/final_bins/ --interleaved {input.pe_1} --min-covered-fraction 0.0 -x fna > bins/coverm_abundances.tsv; "
