rule rerun_rosella:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    group: 'binning'
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


rule benchmark_vamb:
    input:
        "data/vamb_bins/done"
    group: 'binning'
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


rule binner_result:
    input:
        metabat2_done = "data/metabat_bins_2/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        rosella_done = "data/rosella_bins/done",
        vamb_done = "data/vamb_bins/done",
    group: 'binning'
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
        "ln -s ../vamb_bins/bins/*.fna ./ && ls *.fna | parallel 'mv {{}} vamb_bins_{{}}'; rm -f vamb_bins_\*.fna; "
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
    group: 'binning'
    output:
        das_tool_done = "data/das_tool_without_rosella/done"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    shell:
        """
        mkdir -p data/no_rosella/
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/no_rosella/metabat_bins_2.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/no_rosella/metabat_bins_sspec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/no_rosella/metabat_bins_ssens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/no_rosella/metabat_bins_sens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/no_rosella/metabat_bins_spec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/no_rosella/concoct_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/no_rosella/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/no_rosella/vamb_bins.tsv; 
        scaffold2bin_files=$(find data/no_rosella/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g"); 
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} \
         -i $scaffold2bin_files \
         -c {input.fasta} \
         -o data/das_tool_without_rosella/das_tool && \
        touch data/das_tool_without_rosella/done
        """

rule checkm_without_rosella:
    input:
        done = "data/das_tool_without_rosella/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    group: 'binning'
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


rule rosella_benchmark:
    input:
        "data/all_bins/done",
        "data/das_tool_bins/done",
        "data/checkm.out",
        "data/checkm_without_rosella.out",
        # "data/coverm_abundances.tsv",
    group: 'binning'
    output:
        "data/done"
    shell:
        "touch data/done"


rule reset_benchmark:
    log:
        temp('data/reset')
    shell:
        'rm -rf data/rosella_bins/; '
        'rm -rf data/metabat*/; '
        'rm -rf data/*_bins/; '
        'rm -rf data/das_tool_*/; '
        'rm -rf data/checkm*; '
        'rm -rf data/all_bins/; '
        'rm -rf data/*_cov.tsv; '
        'rm -rf data/coverm.cov; '
        'rm -f data/done; '
        'touch data/reset'

rule reset_rosella:
    log:
        temp('data/reset_rosella')
    shell:
        'rm -rf data/rosella_bins/; '
        'rm -rf data/das_tool_bins/; '
        'rm -rf data/checkm; '
        'rm -rf data/checkm.out; '
        'rm -rf data/gtdbtk; '
        'rm -rf data/all_bins/; '
        'rm -rf data/done; '
        'touch data/reset_rosella; '