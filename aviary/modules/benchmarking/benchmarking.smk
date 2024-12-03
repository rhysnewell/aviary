rule rerun_rosella:
    input:
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"]
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
        "rosella recover -r {input.fasta} -C {input.coverage} -t {threads} -o data/rosella_bins "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} --n-neighbors 100 && "
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
        semibin_done = "data/semibin_bins/done",
    output:
         "data/all_bins/checkm.out"
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
        "ln -s ../semibin_bins/output_bins/*.fa ./ && ls *.fa | parallel 'mv {{}} {{.}}.semibin.fna'; "
        "ln -s ../concoct_bins/*.fa ./ && ls *.fa | parallel 'mv {{}} concoct_{{.}}.fna'; "
        "ln -s ../maxbin2_bins/*.fasta ./ && ls *.fasta | parallel 'mv {{}} maxbin2_{{.}}.fna'; "
        "ln -s ../rosella_bins/*.fna ./; "
        "rm -f \*.fna; "
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna --tab_table ./ checkm > checkm.out; "
        "checkm qa -o 2 --tab_table -f checkm.out checkm/lineage.ms checkm/; "
        "touch done && cd ../../"


rule checkm2_all_bins:
    input:
        checkm1_done = "data/all_bins/checkm.out"
    output:
        checkm2_report = "data/checkm2_all_bins/quality_report.tsv"
    threads:
        config["max_threads"]
    params:
        checkm2_db_folder = config["checkm2_db_folder"]
    conda:
        "../../envs/checkm2.yaml"
    shell:
        "export CHECKM2DB={params.checkm2_db_folder}; "
        "checkm2 predict -i data/all_bins/ -x fna -o data/checkm2_all_bins/ -t {threads} --force"

rule checkm2_das_tool:
    input:
        das_tool_wr = "data/das_tool_bins_no_refine/checkm.out",
        das_tool_nr = "data/das_tool_without_rosella/done"
    output:
        checkm2_report_wr = "data/checkm2_das_tool_wr/quality_report.tsv",
        checkm2_report_nr = "data/checkm2_das_tool_nr/quality_report.tsv"
    threads:
        config["max_threads"]
    params:
        checkm2_db_folder = config["checkm2_db_folder"]
    conda:
        "../../envs/checkm2.yaml"
    shell:
        "export CHECKM2DB={params.checkm2_db_folder}; "
        "checkm2 predict -i data/das_tool_bins_no_refine/das_tool_DASTool_bins/ -x fa -o data/checkm2_das_tool_wr/ -t {threads} --force; "
        "checkm2 predict -i data/das_tool_without_rosella/das_tool_DASTool_bins/ -x fa -o data/checkm2_das_tool_nr/ -t {threads} --force"


rule checkm2_refined:
    input:
        m2_refined = "data/rosella_refine_metabat2/checkm.out",
        ro_refined = "data/rosella_refine_rosella/checkm.out",
        sb_refined = "data/rosella_refine_semibin/checkm.out",
        dt_refined = "data/rosella_refine_das_tool/checkm.out"
    output:
        m2_report = "data/m2_refined_checkm2/quality_report.tsv",
        ro_report = "data/ro_refined_checkm2/quality_report.tsv",
        sb_report = "data/sb_refined_checkm2/quality_report.tsv",
        dt_report = "data/dt_refined_checkm2/quality_report.tsv"
    threads:
        config["max_threads"]
    params:
        checkm2_db_folder = config["checkm2_db_folder"]
    conda:
        "../../envs/checkm2.yaml"
    shell:
        "export CHECKM2DB={params.checkm2_db_folder}; "
        "checkm2 predict -i data/rosella_refine_metabat2/final_bins/ -x fna -o data/m2_refined_checkm2/ -t {threads} --force || touch {output.m2_report}; "
        "checkm2 predict -i data/rosella_refine_rosella/final_bins/ -x fna -o data/ro_refined_checkm2/ -t {threads} --force || touch {output.ro_report}; "
        "checkm2 predict -i data/rosella_refine_semibin/final_bins/ -x fna -o data/sb_refined_checkm2/ -t {threads} --force || touch {output.sb_report}; "
        "checkm2 predict -i data/rosella_refine_das_tool/final_bins/ -x fna -o data/dt_refined_checkm2/ -t {threads} --force || touch {output.dt_report}; "

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
        checkm = 'data/rosella_bins/checkm.out',
        rosella = 'data/rosella_bins/done',
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/rosella_refine_rosella/done'
    params:
        bin_folder = "data/rosella_bins/",
        extension = "fna",
        output_folder = "data/rosella_refine_rosella/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 1,
        max_retries = config["refinery_max_retries"],
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 10,
        final_refining = False,
        bin_prefix = "rosella"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/rosella.yaml"
    script:
        "../binning/scripts/rosella_refine.py"

rule rosella_refine_benchmark_2:
    input:
        checkm = 'data/metabat_bins_2/checkm.out',
        rosella = 'data/metabat_bins_2/done',
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/rosella_refine_metabat2/done'
    params:
        bin_folder = "data/metabat_bins_2/",
        extension = "fa",
        output_folder = "data/rosella_refine_metabat2/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 1,
        max_retries = config["refinery_max_retries"],
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 10,
        final_refining = False,
        bin_prefix = "metabat2"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/rosella.yaml"
    script:
        "../binning/scripts/rosella_refine.py"

rule rosella_refine_benchmark_3:
    input:
        checkm = 'data/semibin_bins/checkm.out',
        rosella = 'data/semibin_bins/done',
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/rosella_refine_semibin/done'
    params:
        bin_folder = "data/semibin_bins/output_bins/",
        extension = "fa",
        output_folder = "data/rosella_refine_semibin/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 1,
        max_retries = config["refinery_max_retries"],
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 10,
        final_refining = False,
        bin_prefix = "semibin2"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/rosella.yaml"
    script:
        "../binning/scripts/rosella_refine.py"


rule rosella_refine_benchmark_4:
    input:
        checkm = 'data/das_tool_bins_with_refine/checkm.out',
        das_tool = 'data/das_tool_bins_with_refine/done',
        coverage = "data/coverm.cov",
        fasta = config["fasta"],
        # kmers = "data/rosella_bins/rosella_kmer_table.tsv"
    output:
        'data/rosella_refine_das_tool/done'
    params:
        bin_folder = "data/das_tool_bins_pre_refine/das_tool_DASTool_bins/",
        extension = "fa",
        output_folder = "data/rosella_refine_das_tool/",
        min_bin_size = config["min_bin_size"],
        max_iterations = 1,
        max_retries = config["refinery_max_retries"],
        pplacer_threads = config["pplacer_threads"],
        max_contamination = 10,
        final_refining = True,
        bin_prefix = "dastool"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/rosella.yaml"
    script:
        "../binning/scripts/rosella_refine.py"


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
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna "
        "--tab_table data/rosella_refine_rosella/final_bins/ data/rosella_refine_rosella/checkm > {output.checkm1_done} || touch {output.checkm1_done}; "
        "checkm qa -o 2 --tab_table -f {output.checkm1_done} "
        "data/rosella_refine_rosella/checkm/lineage.ms data/rosella_refine_rosella/checkm/; "

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
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna "
        "--tab_table data/rosella_refine_metabat2/final_bins/ data/rosella_refine_metabat2/checkm > {output.checkm1_done} || touch {output.checkm1_done}; "
        "checkm qa -o 2 --tab_table -f {output.checkm1_done} "
        "data/rosella_refine_metabat2/checkm/lineage.ms data/rosella_refine_metabat2/checkm/; "

rule checkm1_rosella_refine_3:
    input:
        refine1 = "data/rosella_refine_semibin/done",
    output:
        checkm1_done = "data/rosella_refine_semibin/checkm.out",
    params:
        pplacer_threads = config["pplacer_threads"]
    threads:
        config["max_threads"]
    conda:
        "../../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna "
        "--tab_table data/rosella_refine_semibin/final_bins/ data/rosella_refine_semibin/checkm > {output.checkm1_done} || touch {output.checkm1_done}; "
        "checkm qa -o 2 --tab_table -f {output.checkm1_done} "
        "data/rosella_refine_semibin/checkm/lineage.ms data/rosella_refine_semibin/checkm/; "

rule checkm1_rosella_refine_4:
    input:
        refine1 = "data/rosella_refine_das_tool/done",
    output:
        checkm1_done = "data/rosella_refine_das_tool/checkm.out",
    params:
        pplacer_threads = config["pplacer_threads"]
    threads:
        config["max_threads"]
    conda:
        "../../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} -x fna "
        "--tab_table data/rosella_refine_das_tool/ data/rosella_refine_das_tool/checkm > {output.checkm1_done} || touch {output.checkm1_done}; "
        "checkm qa -o 2 --tab_table -f {output.checkm1_done} "
        "data/rosella_refine_das_tool/checkm/lineage.ms data/rosella_refine_das_tool/checkm/; "

rule das_tool_with_refine:
    """
    Runs dasTool on the output of all binning algorithms. If a binner failed to produce bins then their output is ignored
    """
    input:
        fasta = config["fasta"],
        rosella_refine = "data/rosella_refine_rosella/done",
        metabat_refine = "data/rosella_refine_metabat2/done",
        semibin_refine = "data/rosella_refine_semibin/done",
        # metabat2_done = "data/metabat2_refined/done",
        concoct_done = "data/concoct_bins/done",
        maxbin_done = "data/maxbin2_bins/done",
        metabat_sspec = "data/metabat_bins_sspec/done",
        metabat_spec = "data/metabat_bins_spec/done",
        metabat_ssens = "data/metabat_bins_ssens/done",
        metabat_sense = "data/metabat_bins_sens/done",
        # rosella_done = "data/rosella_refined/done",
        vamb_done = "data/vamb_bins/done"
    output:
        das_tool_done = "data/das_tool_bins_with_refine/done"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    benchmark:
        "benchmarks/das_tool_refine.benchmark.txt"
    shell:
        """
        mkdir -p data/refine_input; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_refine_rosella/final_bins/ -e fna > data/refine_input/rosella_refine_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_refine_metabat2/final_bins/ -e fna > data/refine_input/metabat_refine_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_refine_semibin/final_bins/ -e fna > data/refine_input/semibin_refine_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/refine_input/metabat_bins_sspec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/refine_input/metabat_bins_ssens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/refine_input/metabat_bins_sens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/refine_input/metabat_bins_spec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/refine_input/concoct_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/refine_input/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/refine_input/vamb_bins.tsv; 
        scaffold2bin_files=$(find data/refine_input/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g"); 
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} --score_threshold -42 \
         -i $scaffold2bin_files \
         -c {input.fasta} \
         -o data/das_tool_bins_with_refine/das_tool && \
        touch data/das_tool_bins_with_refine/done
        """

rule checkm_das_tool_refine:
    input:
        done = "data/das_tool_bins_with_refine/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/das_tool_bins_with_refine/checkm.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_bins_with_refine/das_tool_DASTool_bins data/das_tool_bins_with_refine/checkm --tab_table -f {output[0]}; '
        'checkm qa -o 2 --tab_table -f {output[0]} '
        'data/das_tool_bins_with_refine/checkm/lineage.ms data/das_tool_bins_with_refine/checkm/; '

rule benchmark_refine:
    input:
        "data/rosella_refine_metabat2/checkm.out",
        "data/rosella_refine_rosella/checkm.out",
        "data/rosella_refine_semibin/checkm.out",
        "data/rosella_refine_das_tool/checkm.out"
    output:
        temp("data/refined")
    shell:
        "touch data/refined"

rule bin_statistics:
    input:
        checkm1_done = "data/all_bins/checkm.out",
        m2_refined = "data/rosella_refine_metabat2/checkm.out",
        ro_refined = "data/rosella_refine_rosella/checkm.out",
        sb_refined = "data/rosella_refine_semibin/checkm.out",
        dt_refined = "data/rosella_refine_das_tool/checkm.out",
        dastool_nr = "data/das_tool_bins_no_refine/checkm.out",
        dastool_wr = "data/das_tool_bins_with_refine/checkm.out",
        dastool_wo = "data/checkm_without_rosella.out",
        coverage_file = "data/coverm.cov"
    output:
        bin_stats = "data/bin_stats/all_bin_stats.tsv",
        m2_stats = "data/bin_stats/m2_refined_bin_stats.tsv",
        ro_stats = "data/bin_stats/ro_refined_bin_stats.tsv",
        sb_stats = "data/bin_stats/sb_refined_bin_stats.tsv",
        dt_stats = "data/bin_stats/dt_refined_bin_stats.tsv",
        das_tool_wr_stats = "data/bin_stats/das_tool_wr_stats.tsv",
        das_tool_nr_stats = "data/bin_stats/das_tool_nr_stats.tsv",
        das_tool_wo_stats = "data/bin_stats/das_tool_wo_stats.tsv",
    run:
        import glob
        import pandas as pd
        from Bio import SeqIO
        import os
        from pathlib import Path

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
        das_tool_wo_checkm = pd.read_csv(input.dastool_wo, sep='\t', comment="[")

        das_tool_wr_stats = retrieve_stats("data/das_tool_bins_with_refine/das_tool_DASTool_bins/", "fa", coverage_file, das_tool_wr_checkm)
        das_tool_wr_stats = pd.DataFrame(data = das_tool_wr_stats)
        das_tool_wr_stats.to_csv(output.das_tool_wr_stats, sep='\t', index=False)

        das_tool_nr_stats = retrieve_stats("data/das_tool_bins_no_refine/das_tool_DASTool_bins/", "fa", coverage_file, das_tool_nr_checkm)
        das_tool_nr_stats = pd.DataFrame(data = das_tool_nr_stats)
        das_tool_nr_stats.to_csv(output.das_tool_nr_stats, sep='\t', index=False)

        das_tool_wo_stats = retrieve_stats("data/das_tool_without_rosella/das_tool_DASTool_bins/", "fa", coverage_file, das_tool_wo_checkm)
        das_tool_wo_stats = pd.DataFrame(data = das_tool_wo_stats)
        das_tool_wo_stats.to_csv(output.das_tool_wo_stats, sep='\t', index=False)

        # repeat for both refine runs
        try:
            ro_checkm = pd.read_csv(input.ro_refined, sep='\t', comment="[")
            ro_stats = retrieve_stats("data/rosella_refine_rosella/", "fna", coverage_file, ro_checkm)
            ro_stats = pd.DataFrame(data = ro_stats)
            ro_stats.to_csv(output.ro_stats, sep='\t', index=False)
        except FileNotFoundError:
            Path(output.ro_stats).touch()

        try:
            m2_checkm = pd.read_csv(input.m2_refined, sep='\t', comment="[")
            m2_stats = retrieve_stats("data/rosella_refine_metabat2/", "fna", coverage_file, m2_checkm)
            m2_stats = pd.DataFrame(data = m2_stats)
            m2_stats.to_csv(output.m2_stats, sep='\t', index=False)
        except FileNotFoundError:
            Path(output.m2_stats).touch()

        try:
            sb_checkm = pd.read_csv(input.sb_refined, sep='\t', comment="[")
            sb_stats = retrieve_stats("data/rosella_refine_semibin/", "fna", coverage_file, sb_checkm)
            sb_stats = pd.DataFrame(data = sb_stats)
            sb_stats.to_csv(output.sb_stats, sep='\t', index=False)
        except FileNotFoundError:
            Path(output.sb_stats).touch()

        try:
            dt_checkm = pd.read_csv(input.dt_refined, sep='\t', comment="[")
            dt_stats = retrieve_stats("data/rosella_refine_das_tool/", "fna", coverage_file, dt_checkm)
            dt_stats = pd.DataFrame(data = dt_stats)
            dt_stats.to_csv(output.dt_stats, sep='\t', index=False)
        except FileNotFoundError:
            Path(output.dt_stats).touch()

rule get_stats:
    input:
        bin_stats = "data/bin_stats/all_bin_stats.tsv",
        m2_refined_stats = "data/bin_stats/m2_refined_bin_stats.tsv",
        ro_refined_stats = "data/bin_stats/ro_refined_bin_stats.tsv",
        sb_refined_stats = "data/bin_stats/sb_refined_bin_stats.tsv",
        dt_refined_stats = "data/bin_stats/dt_refined_bin_stats.tsv",
        das_tool_wr_stats = "data/bin_stats/das_tool_wr_stats.tsv",
        das_tool_wo_stats = "data/bin_stats/das_tool_wo_stats.tsv",
        das_tool_nr_stats = "data/bin_stats/das_tool_nr_stats.tsv",
        checkm2_report = "data/checkm2_all_bins/quality_report.tsv",
        checkm2_report_wr = "data/checkm2_das_tool_wr/quality_report.tsv",
        checkm2_report_nr = "data/checkm2_das_tool_nr/quality_report.tsv",
        m2_report = "data/m2_refined_checkm2/quality_report.tsv",
        ro_report = "data/ro_refined_checkm2/quality_report.tsv",
        sb_report = "data/sb_refined_checkm2/quality_report.tsv",
        dt_report = "data/dt_refined_checkm2/quality_report.tsv"
    output:
        done = "data/stats_done"
    shell:
        "touch data/stats_done"

rule das_tool_without_rosella:
    input:
        fasta = config["fasta"],
        metabat2_done = "data/metabat_bins_2/done",
        semibin_done = "data/semibin_bins/done",
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
    benchmark:
        "benchmarks/das_tool_no_rosella.benchmark.txt"
    conda:
        "../binning/envs/das_tool.yaml"
    shell:
        """
        mkdir -p data/no_rosella/
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2 -e fa > data/no_rosella/metabat_bins_2.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/semibin_bins/output_bins/ -e fa > data/no_rosella/semibin_bins.tsv; 
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

rule das_tool_no_refine:
    """
    Runs dasTool on the output of all binning algorithms. If a binner failed to produce bins then their output is ignored
    """
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
        semibin_done = "data/semibin_bins/done",
        vamb_done = "data/vamb_bins/done",
    output:
        das_tool_done = "data/das_tool_bins_no_refine/done"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    benchmark:
        "benchmarks/das_tool_no_refine.benchmark.txt"
    shell:
        """
        mkdir -p data/no_refine_input/; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/no_refine_input/metabat_bins_sspec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/no_refine_input/metabat_bins_ssens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/no_refine_input/metabat_bins_sens.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/no_refine_input/metabat_bins_spec.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/no_refine_input/concoct_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/no_refine_input/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/no_refine_input/vamb_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/rosella_bins/ -e fna > data/no_refine_input/rosella_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/semibin_bins/output_bins/ -e fa > data/no_refine_input/semibin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2/ -e fa > data/no_refine_input/metabat2_bins.tsv; 
        scaffold2bin_files=$(find data/no_refine_input/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g"); 
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} --score_threshold -42 \
         -i $scaffold2bin_files \
         -c {input.fasta} \
         -o data/das_tool_bins_no_refine/das_tool && \
        touch data/das_tool_bins_no_refine/done
        """

rule checkm_das_tool_no_refine:
    input:
        done = "data/das_tool_bins_no_refine/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/das_tool_bins_no_refine/checkm.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/das_tool_bins_no_refine/das_tool_DASTool_bins data/das_tool_bins_no_refine/checkm --tab_table -f {output[0]}; '
        'checkm qa -o 2 --tab_table -f {output[0]} '
        'data/das_tool_bins_no_refine/checkm/lineage.ms data/das_tool_bins_no_refine/checkm/; '

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
        '--tab_table -f {output[0]}; '
        'checkm qa -o 2 --tab_table -f {output[0]} '
        'data/checkm_without_rosella/lineage.ms data/checkm_without_rosella/; '

rule verify_gsa_mapping_file:
    input:
        gsa_mapping = config["gsa_mappings"]
    output:
        gsa_output = "data/gsa_mapping.binning"
    run:
        good_input = False
        with open(output.gsa_output, 'w') as verified_output:
            with open(input.gsa_mapping, 'r') as gsa_mapping:
                for idx, line_original in enumerate(gsa_mapping):
                    line_original = line_original.strip()
                    if idx == 0 and line_original.startswith("@Version"):
                        # input is already good to go
                        good_input = True
                    if good_input:
                        print(line_original, file=verified_output)
                    else:
                        if idx == 0:
                            # print gsa header
                            print("@Version:0.9.1", file=verified_output)
                            print("@SampleID:gsa", file=verified_output)
                            print("", file=verified_output)
                            print("@@SEQUENCEID\tBINID\tTAXID", file=verified_output)
                        else:
                            line_original = line_original.split('\t')
                            print("\t".join(line_original[0:3]), file=verified_output)

rule add_lengths:
    input:
        # gsa_mappings = config["gsa_mappings"],
        gsa_mappings = "data/gsa_mapping.binning",
        gsa = config["fasta"]
    output:
        "data/gsa_mapping_lengths.binning"
    shell:
        "add_length_column.py -g {input.gsa_mappings} -f {input.gsa} > {output[0]}"


rule amber_for_refining:
    input:
        semibin = "data/semibin_bins/done",
        rosella = "data/rosella_bins/done",
        m2 = "data/metabat_bins_2/done",
        gsa_mappings = "data/gsa_mapping_lengths.binning"
    output:
        amber_refine = "data/amber_refine/for_refine/index.html"
    shell:
        'mkdir -p data/amber_refine; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_refine/rosella_amber.tsv data/rosella_bins/*.fna || rm data/amber_refine/rosella_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_refine/m2_amber.tsv data/metabat_bins_2/*.fa || rm data/amber_refine/m2_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_refine/semibin_amber.tsv data/semibin_bins/output_bins/*.fa || rm data/amber_refine/semibin_amber.tsv; '
        'sed -i "s/_SAMPLEID_/gsa/" data/amber_refine/*.tsv; '
        'amber.py -g {input.gsa_mappings} -o data/amber_refine/for_refine -x 90,80,70,60,50 -y 10,5 '
        'data/amber_refine/*.tsv'


rule setup_amber:
    input:
        m1_sspec = "data/metabat_bins_sspec/done",
        m1_ssens = "data/metabat_bins_ssens/done",
        m1_sens = "data/metabat_bins_sens/done",
        m1_spec = "data/metabat_bins_spec/done",
        concoct = "data/concoct_bins/done",
        maxbin = "data/maxbin2_bins/done",
        vamb = "data/vamb_bins/done",
        semibin = "data/semibin_bins/done",
        rosella = "data/rosella_bins/done",
        m2 = "data/metabat_bins_2/done",
        semibin_refined = "data/rosella_refine_semibin/done",
        rosella_refined = "data/rosella_refine_rosella/done",
        m2_refined = "data/rosella_refine_metabat2/done",
        dt_wor = "data/das_tool_without_rosella/done",
        dt_nr = "data/das_tool_bins_no_refine/done",
        dt_wr = "data/das_tool_bins_with_refine/done",
        gsa_mappings = "data/gsa_mapping_lengths.binning"
        # gsa_mappings = config["gsa_mappings"]
    output:
        amber_ready = "data/amber_out/done"
    shell:
        'mkdir -p data/amber_out; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/rosella_amber.tsv data/rosella_bins/*.fna || rm data/amber_out/rosella_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/m1_sspec_amber.tsv data/metabat_bins_sspec/*.fa || rm data/amber_out/m1_sspec_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/m1_spec_amber.tsv data/metabat_bins_spec/*.fa || rm data/amber_out/m1_spec_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/m1_ssens_amber.tsv data/metabat_bins_ssens/*.fa || rm data/amber_out/m1_ssens_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/m1_sens_amber.tsv data/metabat_bins_sens/*.fa || rm data/amber_out/m1_sens_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/m2_amber.tsv data/metabat_bins_2/*.fa || rm data/amber_out/m2_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/concoct_amber.tsv data/concoct_bins/*.fa || rm data/amber_out/concoct_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/maxbin2_amber.tsv data/maxbin2_bins/*.fasta || rm data/amber_out/maxbin2_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/vamb_amber.tsv data/vamb_bins/bins/*.fna || rm data/amber_out/vamb_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/semibin_amber.tsv data/semibin_bins/output_bins/*.fa || rm data/amber_out/semibin_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/semibin_refined_amber.tsv data/rosella_refine_semibin/final_bins/*.fna || rm data/amber_out/semibin_refined_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/rosella_refined_amber.tsv data/rosella_refine_rosella/final_bins/*.fna || rm data/amber_out/rosella_refined_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/metabat2_refined_amber.tsv data/rosella_refine_metabat2/final_bins/*.fna || rm data/amber_out/metabat2_refined_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/dt_wor_amber.tsv data/das_tool_without_rosella/das_tool_DASTool_bins/*.fa || rm data/amber_out/dt_wor_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/dt_nr_amber.tsv data/das_tool_bins_no_refine/das_tool_DASTool_bins/*.fa || rm data/amber_out/dt_nr_amber.tsv; '
        'convert_fasta_bins_to_biobox_format.py -o data/amber_out/dt_wr_amber.tsv data/das_tool_bins_with_refine/das_tool_DASTool_bins/*.fa || rm data/amber_out/dt_wr_amber.tsv; '
        'touch data/amber_out/done; '

rule run_amber:
    input:
        amber_ready = "data/amber_out/done",
        gsa_mappings = "data/gsa_mapping_lengths.binning"
    output:
        amber_out = "data/amber_out/benchmarks/index.html"
    # "scaffold2bin_files=$(find data/amber_out/*.tsv -not -empty -exec ls {{}} \; | tr '\n' ',' | sed 's/,$//g'); "
    shell:
        "sed -i 's/_SAMPLEID_/gsa/' data/amber_out/*.tsv; "
        "amber.py -g {input.gsa_mappings} -o data/amber_out/benchmarks -x 90,80,70,60,50 -y 10,5 "
        "data/amber_out/*.tsv"

rule rosella_benchmark:
    input:
        "data/all_bins/checkm.out",
        "data/das_tool_bins_no_refine/checkm.out",
        "data/checkm_without_rosella.out",
        # "data/coverm_abundances.tsv",
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
        'rm -rf data/rosella_refine*/; '
        'rm -rf data/no_rosella/; '
        'rm -rf data/das_tool_*/; '
        'rm -rf data/checkm_*; '
        'rm -rf data/checkm.out; '
        'rm -rf data/gtdbtk; '
        'rm -rf data/all_bins/; '
        'rm -rf data/done; '
        'rm -rf bins/; '
        'rm -rf taxonomy/; '
        'rm -rf diversity/'
        'touch data/reset_rosella; '

rule reset_recover:
    log:
        temp('data/reset_recover')
    shell:
        'rm -f data/rosella_bins.tsv; '
        'rm -f data/metabat_bin_2.tsv; '
        'rm -rf data/refined_bins/; '
        'rm -rf data/rosella_refined/; '
        'rm -rf data/metabat2_refined/; '
        'rm -rf data/das_tool_bins/; '
        'rm -rf data/checkm; '
        'rm -rf data/checkm.out; '
        'rm -rf data/gtdbtk; '
        'rm -rf data/all_bins/; '
        'rm -rf data/done; '
        'rm -rf bins/; '
        'rm -rf taxonomy/; '
        'rm -rf diversity/'
        'touch data/reset_recover; '

rule magpurify_metabat2:
    input:
        "data/metabat_bins_2/done",
    output:
        output_folder = directory("data/magpurify_metabat2/"),
        checkm_out = "data/magpurify_metabat2/cleaned_bins/checkm.out"
    threads:
        config["max_threads"]
    conda:
        "envs/magpurify.yaml"
    params:
        pplacer_threads = config["pplacer_threads"]
    shell:
        "export MAGPURIFYDB=/home/n10853499/databases/MAGpurify-db-v1.0; "
        "mkdir -p {output.output_folder}; "
        "mkdir -p {output.output_folder}/cleaned_bins/; "
        "ls data/metabat_bins_2/*.fa | parallel -j {threads} "
        "'magpurify phylo-markers {{}} {output.output_folder}/{{/.}}; "
        "magpurify clade-markers {{}} {output.output_folder}/{{/.}}; "
        "magpurify tetra-freq {{}} {output.output_folder}/{{/.}}; "
        "magpurify gc-content {{}} {output.output_folder}/{{/.}}; "
        "magpurify known-contam {{}} {output.output_folder}/{{/.}}; "
        "magpurify clean-bin {{}} {output.output_folder}/{{/.}} {output.output_folder}/cleaned_bins/{{/.}}_cleaned.fna'; "
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} "
        "-x fna {output.output_folder}/cleaned_bins/ {output.output_folder}/cleaned_bins/checkm "
        "--tab_table -f {output.output_folder}/cleaned_bins/checkm.out; "
        "checkm qa -o 2 --tab_table -f {output.output_folder}/cleaned_bins/checkm.out "
        "{output.output_folder}/cleaned_bins/checkm/lineage.ms {output.output_folder}/cleaned_bins/checkm/; "

rule magpurify_rosella:
    input:
        "data/rosella_bins/done",
    output:
        output_folder = directory("data/magpurify_rosella/"),
        checkm_out = "data/magpurify_rosella/cleaned_bins/checkm.out"
    threads:
        config["max_threads"]
    conda:
        "envs/magpurify.yaml"
    params:
        pplacer_threads = config["pplacer_threads"]
    shell:
        "export MAGPURIFYDB=/home/n10853499/databases/MAGpurify-db-v1.0; "
        "mkdir -p {output.output_folder}; "
        "mkdir -p {output.output_folder}/cleaned_bins/; "
        "ls data/rosella_bins/*.fna | parallel -j {threads} "
        "'magpurify phylo-markers {{}} {output.output_folder}/{{/.}}; "
        "magpurify clade-markers {{}} {output.output_folder}/{{/.}}; "
        "magpurify tetra-freq {{}} {output.output_folder}/{{/.}}; "
        "magpurify gc-content {{}} {output.output_folder}/{{/.}}; "
        "magpurify known-contam {{}} {output.output_folder}/{{/.}}; "
        "magpurify clean-bin {{}} {output.output_folder}/{{/.}} {output.output_folder}/cleaned_bins/{{/.}}_cleaned.fna'; "
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} "
        "-x fna {output.output_folder}/cleaned_bins/ {output.output_folder}/cleaned_bins/checkm "
        "--tab_table -f {output.output_folder}/cleaned_bins/checkm.out; "
        "checkm qa -o 2 --tab_table -f {output.output_folder}/cleaned_bins/checkm.out "
        "{output.output_folder}/cleaned_bins/checkm/lineage.ms {output.output_folder}/cleaned_bins/checkm/; "

rule complete_magpurify:
    input:
         checkm_out_1 = "data/magpurify_metabat2/cleaned_bins/checkm.out",
         checkm_out_2 = "data/magpurify_rosella/cleaned_bins/checkm.out"
    log:
       temp("data/magpurify_done")

rule atlas_dastool:
    input:
         fasta = config["fasta"],
         m2_done = "data/metabat_bins_2/done",
         maxbin_done = "data/maxbin2_bins/done"
    output:
        das_tool_done = "data/atlas_dastool/done"
    threads:
        config["max_threads"]
    conda:
        "../binning/envs/das_tool.yaml"
    benchmark:
        "benchmarks/das_tool.benchmark.txt"
    shell:
        """
        Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/maxbin_bins.tsv; 
        Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2/ -e fa > data/metabat2_bins.tsv;
        DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} --score_threshold -42 \
         -i data/maxbin_bins.tsv,data/metabat2_bins.tsv \
         -c {input.fasta} \
         -o data/atlas_dastool/das_tool && \
        touch data/atlas_dastool/done
        """

rule simulate_atlas:
    input:
         dastool_done = "data/atlas_dastool/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/atlas_dastool/checkm.out"
    conda:
        "../../envs/checkm.yaml"
    threads:
        config["max_threads"]
    shell:
        'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
        '-x fa data/atlas_dastool/das_tool_DASTool_bins data/atlas_dastool/checkm --tab_table -f {output[0]}; '
        'checkm qa -o 2 --tab_table -f {output[0]} '
        'data/atlas_dastool/checkm/lineage.ms data/atlas_dastool/checkm/; '


# rule dastool_no_refine:
#     """
#     Runs dasTool on the output of all binning algorithms. If a binner failed to produce bins then their output is ignored
#     Does not include refined results
#     """
#     input:
#         fasta = config["fasta"],
#         metabat2_done = "data/metabat_bins_2/done",
#         concoct_done = "data/concoct_bins/done",
#         maxbin_done = "data/maxbin2_bins/done",
#         metabat_sspec = "data/metabat_bins_sspec/done",
#         metabat_spec = "data/metabat_bins_spec/done",
#         metabat_ssens = "data/metabat_bins_ssens/done",
#         metabat_sense = "data/metabat_bins_sens/done",
#         rosella_done = "data/rosella_bins/done",
#         vamb_done = "data/vamb_bins/done",
#     output:
#         das_tool_done = "data/das_tool_bins_no_refine/done"
#     threads:
#         config["max_threads"]
#     conda:
#         "../binning/envs/das_tool.yaml"
#     benchmark:
#         "benchmarks/das_tool.benchmark.txt"
#     shell:
#         """
#         mkdir -p data/unrefined;
#         Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sspec -e fa > data/unrefined/metabat_bins_sspec.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_ssens -e fa > data/unrefined/metabat_bins_ssens.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_sens -e fa > data/unrefined/metabat_bins_sens.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_spec -e fa > data/unrefined/metabat_bins_spec.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/concoct_bins -e fa > data/unrefined/concoct_bins.tsv || touch data/unrefined/concoct_bins.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/maxbin2_bins -e fasta > data/unrefined/maxbin_bins.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/vamb_bins/bins -e fna > data/unrefined/vamb_bins.tsv || touch data/unrefined/vamb_bins.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/rosella_bins/ -e fna > data/unrefined/rosella_bins.tsv;
#         Fasta_to_Scaffolds2Bin.sh -i data/metabat_bins_2/ -e fa > data/unrefined/metabat2_bins.tsv;
#         scaffold2bin_files=$(find data/unrefined/*bins*.tsv -not -empty -exec ls {{}} \; | tr "\n" ',' | sed "s/,$//g");
#         DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {threads} --score_threshold -42 \
#          -i $scaffold2bin_files \
#          -c {input.fasta} \
#          -o data/das_tool_bins_no_refine/das_tool && \
#         touch data/das_tool_bins_no_refine/done
#         """

# rule checkm_no_refine:
#     input:
#          dastool_done = "data/das_tool_bins_no_refine/done"
#     params:
#         pplacer_threads = config["pplacer_threads"]
#     output:
#         "data/das_tool_bins_no_refine/checkm.out"
#     conda:
#         "../../envs/checkm.yaml"
#     threads:
#         config["max_threads"]
#     shell:
#         'checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} '
#         '-x fa data/das_tool_bins_no_refine/das_tool_DASTool_bins data/das_tool_bins_no_refine/checkm --tab_table -f data/das_tool_bins_no_refine/checkm.out'
