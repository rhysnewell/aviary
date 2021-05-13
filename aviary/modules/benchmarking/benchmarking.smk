rule rerun_rosella:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    output:
        "data/rosella_bins/rerun"
    conda:
        "../binning/envs/rosella.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/rosella_rerun.benchmark.txt"
    shell:
        "rm -f data/rosella_bins/*.fna; rm -f data/rosella_bins/checkm.out; rm -rf data/rosella_bins/checkm/; "
        "rosella bin -r {input.fasta} -i {input.coverage} -t {threads} -o data/rosella_bins --b-tail 0.4 "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} && "
        "touch data/rosella_bins/rerun"


rule rosella_checkm:
    input:
        done = "data/rosella_bins/rerun"
    params:
        pplacer_threads = config['pplacer_threads']
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fna data/rosella_bins/ data/rosella_bins/checkm '
        '--tab_table -f data/rosella_bins/checkm.out && rm data/rosella_bins/rerun'


rule das_tool_without_rosella:
    input:
        fasta = config["fasta"],
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        vamb_done = "data/vamb_bins/done"
    output:
        das_tool_done = "data/das_tool_without_rosella/done"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    shell:
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/metabat_bins_2.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/metabat_bins_spec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins/ -e fna > data/vamb_bins.tsv && "
        "DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads}"
        " -i data/metabat_bins_2.tsv,data/vamb_bins.tsv,data/metabat_bins_sspec.tsv,data/metabat_bins_spec.tsv,data/metabat_bins_ssens.tsv,data/metabat_bins_sens.tsv,data/maxbin_bins.tsv,data/concoct_bins.tsv"
        " -c {input.fasta} -o data/das_tool_without_rosella/das_tool && "
        "touch data/das_tool_without_rosella/done"

rule das_tool_without_rosella_no_vamb:
    input:
        fasta = config["fasta"],
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        skipped_vamb = "data/vamb_bins/skipped"
    output:
        vamb_skipped = "data/das_tool_without_rosella/skipped_vamb"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    shell:
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/metabat_bins_2.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/metabat_bins_sspec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/metabat_bins_ssens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/metabat_bins_sens.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/metabat_bins_spec.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/concoct_bins.tsv && "
        "Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv && "
        "DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads}"
        " -i data/metabat_bins_2.tsv,data/metabat_bins_sspec.tsv,data/metabat_bins_spec.tsv,data/metabat_bins_ssens.tsv,data/metabat_bins_sens.tsv,data/maxbin_bins.tsv,data/concoct_bins.tsv"
        " -c {input.fasta} -o data/das_tool_without_rosella/das_tool && "
        "touch data/das_tool_without_rosella/done; "
        "touch data/das_tool_without_rosella/skipped_vamb"


rule checkm_without_rosella:
    input:
        done = "data/das_tool_without_rosella/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/checkm_without_rosella.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_without_rosella/das_tool_DASTool_bins data/checkm_without_rosella '
        '--tab_table -f data/checkm_without_rosella.out; '


rule checkm_without_rosella_no_vamb:
    input:
        done = "data/das_tool_without_rosella/skipped_vamb"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/checkm_without_rosella/skipped_vamb"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_without_rosella/das_tool_DASTool_bins data/checkm '
        '--tab_table -f data/checkm_without_rosella.out; '
        'touch data/checkm_without_rosella/skipped_vamb'


rule rosella_benchmark:
    input:
        "data/all_bins/done",
        "data/das_tool_bins/done",
        "data/checkm.out",
        "data/checkm_without_rosella.out",
        # "data/coverm_abundances.tsv",
    output:
        "data/done"
    shell:
        "touch data/done"


rule rosella_benchmark_no_vamb:
    input:
        "data/all_bins/skipped_vamb",
        "data/vamb_bins/skipped",
        "data/checkm/skipped_vamb",
        "data/checkm_without_rosella/skipped_vamb",
        "data/das_tool_without_rosella/skipped_vamb",
        "data/das_tool_bins/skipped_vamb"
    output:
        "data/skipped_vamb"
    shell:
        "touch data/skipped_vamb"

rule reset_benchmark:
    log:
        temp('data/reset')
    shell:
        'rm -rf data/rosella_bins/; '
        'rm -rf data/metabat_bins_2/; '
        'rm -rf data/das_tool_*/; '
        'rm -rf data/checkm*; '
        'rm -rf data/all_bins/; '
        'touch data/reset_all'

rule reset_rosella:
    log:
        temp('data/reset_rosella')
    shell:
        'rm -rf data/rosella_bins/; '
        'rm -rf data/das_tool_bins/; '
        'rm -rf data/checkm; '
        'rm -rf data/checkm.out'
        'rm -rf data/all_bins/; '
        'touch data/reset_rosella; '