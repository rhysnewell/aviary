rule rerun_rosella:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
    group: 'binning'
    output:
        temp("data/rosella_bins/rerun")
    conda:
        "../binning/envs/rosella.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/rosella_rerun.benchmark.txt"
    shell:
        "rm -f data/rosella_bins/*.fna; rm -f data/rosella_bins/checkm.out; rm -rf data/rosella_bins/checkm/; "
        "rosella bin -r {input.fasta} -i {input.coverage} -t {threads} -o data/rosella_bins "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} --n-neighbors 200 && "
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
        '--tab_table -f data/rosella_bins/checkm.out'


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


rule checkm2_all_bins:
    input:
        checkm1_done = "data/all_bins/done"
    output:
        checkm2_report = "data/checkm2_all_bins/quality_report.tsv"
    threads:
        config["max_threads"]
    params:
        checkm2_db_path = config["checkm2_db_path"]
    conda:
        "envs/checkm2.yaml"
    shell:
        "export CHECKM2DB={params.checkm2_db_path}; "
        "checkm2 predict -i data/all_bins/ -x fna -o data/checkm2_all_bins/ -t {threads} --force"

rule checkm2_das_tool:
    input:
        das_tool_wr = "data/das_tool_bins/done",
        das_tool_nr = "data/das_tool_without_rosella/done"
    output:
        checkm2_report_wr = "data/checkm2_das_tool_wr/quality_report.tsv",
        checkm2_report_nr = "data/checkm2_das_tool_nr/quality_report.tsv"
    threads:
        config["max_threads"]
    params:
        checkm2_db_path = config["checkm2_db_path"]
    conda:
        "envs/checkm2.yaml"
    shell:
        "export CHECKM2DB={params.checkm2_db_path}; "
        "checkm2 predict -i data/das_tool_bins/das_tool_DASTool_bins/ -x fa -o data/checkm2_das_tool_wr/ -t {threads} --force; "
        "checkm2 predict -i data/das_tool_without_rosella/das_tool_DASTool_bins/ -x fa -o data/checkm2_das_tool_nr/ -t {threads} --force"

# rule fraction_recovered:
#     input:
#         metabat2_done = "data/metabat_bins_2/done",
#         concoct_done = "data/concoct_bins/done",
#         maxbin_done = "data/maxbin2_bins/done",
#         metabat_sspec = "data/metabat_bins_sspec/done",
#         metabat_spec = "data/metabat_bins_spec/done",
#         metabat_ssens = "data/metabat_bins_ssens/done",
#         metabat_sense = "data/metabat_bins_sens/done",
#         rosella_done = "data/rosella_bins/done",
#         vamb_done = "data/vamb_bins/done"
#     output:
#         metabat2_cov = "data/metabat_bins_2/coverm.cov",
#         concoct_cov = "data/concoct_bins/coverm.cov",
#         maxbin_cov = "data/maxbin2_bins/coverm.cov",
#         metabat_sspec_cov = "data/metabat_bins_sspec/coverm.cov",
#         metabat_spec_cov = "data/metabat_bins_spec/coverm.cov",
#         metabat_ssens_cov = "data/metabat_bins_ssens/coverm.cov",
#         metabat_sense_cov = "data/metabat_bins_sens/coverm.cov",
#         rosella_cov = "data/rosella_bins/coverm.cov",
#         vamb_cov = "data/vamb_bins/coverm.cov"
#     conda:
#         "../../envs/coverm.yaml"
#     shell:
#         ""

rule rosella_refine_benchmark_1:
    input:
        checkm1_done = "data/all_bins/checkm.out",
        coverages = "data/coverm.cov",
        rosella_done = "data/rosella_bins/done",
        assembly = config["fasta"]
    output:
        refine1 = "data/rosella_refine_rosella/done",
    conda:
        "../binning/envs/rosella.yaml",
    threads:
        config["max_threads"]
    shell:
        "rosella refine -a {input.assembly} --coverage-values {input.coverages} -f data/all_bins/rosella*.fna "
        "--checkm-file {input.checkm1_done} -o data/rosella_refine_rosella -t {threads} --contaminated-only; "
        "touch data/rosella_refine_rosella/done; "

rule rosella_refine_benchmark_2:
    input:
        checkm1_done = "data/all_bins/checkm.out",
        coverages = "data/coverm.cov",
        rosella_done = "data/metabat_bins_2/done",
        assembly = config["fasta"]
    output:
        refine2 = "data/rosella_refine_metabat2/done",
    conda:
        "../binning/envs/rosella.yaml",
    threads:
        config["max_threads"]
    shell:
        "rosella refine -a {input.assembly} --coverage-values {input.coverages} -f data/all_bins/*metabat2*.fna "
        "--checkm-file {input.checkm1_done} -o data/rosella_refine_metabat2 -t {threads} --contaminated-only; "
        "touch {output.refine2}; "

rule checkm1_rosella_refine_1:
    input:
        refine1 = "data/rosella_refine_rosella/done",
    output:
        checkm1_done = "data/rosella_refine_rosella/checkm.out",
    params:
        pplacer_threads = config["pplacer_threads"]
    threads:
        config["max_threads"]
    conda:
        "../../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna --tab_table data/rosella_refine_rosella/ data/rosella_refine_rosella/checkm > {output.checkm1_done};"


rule checkm1_rosella_refine_2:
    input:
        refine1 = "data/rosella_refine_metabat2/done",
    output:
        checkm1_done = "data/rosella_refine_metabat2/checkm.out",
    params:
        pplacer_threads = config["pplacer_threads"]
    threads:
        config["max_threads"]
    conda:
        "../../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna --tab_table data/rosella_refine_metabat2/ data/rosella_refine_metabat2/checkm > {output.checkm1_done};"


rule benchmark_refine:
    input:
        "data/rosella_refine_metabat2/checkm.out",
        "data/rosella_refine_rosella/checkm.out"
    output:
        temp("data/refined")
    shell:
        "touch data/refined"

rule bin_statistics:
    input:
        checkm1_done = "data/all_bins/checkm.out",
        dastool_wr = "data/checkm.out",
        dastool_nr = "data/checkm_without_rosella.out",
        coverage_file = "data/coverm.cov"
    output:
        bin_stats = "data/bin_stats/all_bin_stats.tsv",
        das_tool_wr_stats = "data/bin_stats/das_tool_wr_stats.tsv",
        das_tool_nr_stats = "data/bin_stats/das_tool_nr_stats.tsv",
    run:
        import glob
        import pandas as pd
        from Bio import SeqIO
        import os

        if not os.path.exists("data/bin_stats"): os.makedirs("data/bin_stats")

        def retrieve_stats(bin_folder, bin_extension, coverage_file, checkm_file):
            """
            Retrieves information about all bins within `bin_folder` given the file extension `bin_extension`
            returns a dict containing the bin_id, number of contigs, size in base pairs, completeness, and contamination
            """
            output_dict = {"bin_id": [], "n_contigs": [], "size": [], "completeness": [], "contamination": []}

            for fasta_path in glob.glob(f"{bin_folder}/*.{bin_extension}"):
                bin_id = fasta_path.split("/")[-1]
                bin_id = bin_id.replace(f".{bin_extension}", "")
                output_dict["bin_id"].append(bin_id)

                contig_ids = []
                for sequence in SeqIO.parse(open(fasta_path), "fasta"):
                    contig_ids.append(sequence.id)

                output_dict["n_contigs"].append(len(contig_ids))
                output_dict["size"].append(coverage_file[coverage_file["contigName"].isin(contig_ids)]['contigLen'].sum())

                try:
                    # checkm1 uses Bin Id
                    checkm_stats = checkm_file[checkm_file["Bin Id"] == bin_id]
                except KeyError:
                    # checkm2 uses Name
                    checkm_stats = checkm_file[checkm_file["Name"] == bin_id]

                output_dict["completeness"].append(checkm_stats["Completeness"].values[0])
                output_dict["contamination"].append(checkm_stats["Contamination"].values[0])

            return output_dict

        coverage_file = pd.read_csv(input.coverage_file, sep='\t')

        # checkm file for all bins
        all_bins_checkm = pd.read_csv(input.checkm1_done, sep='\t', comment="[")

        # retrieve stats for all bins
        all_bin_stats = retrieve_stats("data/all_bins/", "fna", coverage_file, all_bins_checkm)
        all_bin_stats = pd.DataFrame(data = all_bin_stats)
        all_bin_stats.to_csv(output.bin_stats, sep='\t', index=False)

        # repeat for both das_tool runs
        das_tool_wr_checkm = pd.read_csv(input.dastool_wr, sep='\t', comment="[")
        das_tool_nr_checkm = pd.read_csv(input.dastool_nr, sep='\t', comment="[")

        das_tool_wr_stats = retrieve_stats("data/das_tool_bins/das_tool_DASTool_bins/", "fa", coverage_file, das_tool_wr_checkm)
        das_tool_wr_stats = pd.DataFrame(data = das_tool_wr_stats)
        das_tool_wr_stats.to_csv(output.das_tool_wr_stats, sep='\t', index=False)

        das_tool_nr_stats = retrieve_stats("data/das_tool_without_rosella/das_tool_DASTool_bins/", "fa", coverage_file, das_tool_nr_checkm)
        das_tool_nr_stats = pd.DataFrame(data = das_tool_nr_stats)
        das_tool_nr_stats.to_csv(output.das_tool_nr_stats, sep='\t', index=False)

rule get_stats:
    input:
        bin_stats = "data/bin_stats/all_bin_stats.tsv",
        das_tool_wr_stats = "data/bin_stats/das_tool_wr_stats.tsv",
        das_tool_nr_stats = "data/bin_stats/das_tool_nr_stats.tsv",
        checkm2_report = "data/checkm2_all_bins/quality_report.tsv",
        checkm2_report_wr = "data/checkm2_das_tool_wr/quality_report.tsv",
        checkm2_report_nr = "data/checkm2_das_tool_nr/quality_report.tsv"
    output:
        done = "stats_done"
    shell:
        "touch data/stats_done"

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