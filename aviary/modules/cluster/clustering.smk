onstart:
    import os
    import sys

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


    if config["previous_runs"] == "none":
        sys.exit("Need a list of previous runs to combine for dereplication process")

import os

rule generate_combined_checkm_output:
    output:
        combined_checkm_file = "data/combined_checkm.tsv",
        combined_checkm_stats = "data/combined_checkm_full.tsv",
    params:
        previous_runs = list(config["previous_runs"]),
        use_checkm2_scores = config["use_checkm2_scores"]
    run:
        import pandas as pd

        previous_runs = params.previous_runs
        runs = [pd.read_csv(os.path.abspath(run) + "/bins/checkm_minimal.tsv", sep="\t") for run in previous_runs]
        infos = [pd.read_csv(os.path.abspath(run) + "/bins/bin_info.tsv", sep="\t") for run in previous_runs]

        for (df, run, info_df) in zip(runs, previous_runs, infos):
            df['Bin Id'] = df['Bin Id'].apply(lambda bin: os.path.abspath(run) + "/bins/final_bins/" + str(bin))
            info_df['Bin Id'] = info_df['Bin Id'].apply(lambda bin: os.path.abspath(run) + "/bins/final_bins/" + str(bin))

        concat = pd.concat(runs)
        # need opposite column names to the ones we actually use here
        score_cols = ["Completeness (CheckM2)", "Contamination (CheckM2)"]
        rename_cols = ["Completeness (CheckM1)", "Contamination (CheckM1)"]
        if params.use_checkm2_scores:
            score_cols = ["Completeness (CheckM1)", "Contamination (CheckM1)"]
            rename_cols = ["Completeness (CheckM2)", "Contamination (CheckM2)"]

        concat = concat.loc[:, [col not in score_cols for col in concat.columns]]
        concat.rename({rename_cols[0]: "Completeness", rename_cols[1]: "Contamination"}, inplace=True, axis=1)
        concat.to_csv("data/combined_checkm.tsv", sep="\t", index=False)

        concat_infos = pd.concat(infos)
        concat_infos.to_csv("data/combined_checkm_full.tsv", sep="\t", index=False)


rule generate_genome_list:
    output:
        genome_list = "data/genome_paths.txt"
    params:
        previous_runs = config["previous_runs"]
    shell:
        "for f in {params.previous_runs}; "
        "do "
        "  find $f/bins/final_bins/*.fna; "
        "done >> data/genome_paths.txt"


rule run_galah:
    input:
        genome_list = "data/genome_paths.txt",
        checkm = "data/combined_checkm.tsv"
    output:
        dereplicated_set = "data/dereplicated_clusters.txt"
    params:
        derep_ani = config["ani"],
        precluster_ani = config["precluster_ani"],
        precluster_method = config["precluster_method"],
        min_completeness = config["min_completeness"],
        max_contamination = config["max_contamination"],
    threads:
        config['max_threads']
    conda:
        "../../envs/coverm.yaml"
    shell:
        "galah cluster -t {threads} --checkm-tab-table {input.checkm} " 
        "--genome-fasta-list {input.genome_list} --output-cluster-definition {output.dereplicated_set} "
        "--ani {params.derep_ani} --precluster-ani {params.precluster_ani} --precluster-method {params.precluster_method} "
        "{params.min_completeness} {params.max_contamination} "
        "--output-representative-fasta-directory representative_genomes/ "

rule representative_checkm:
    input:
        dereplicated_clusters = "data/dereplicated_clusters.txt",
        combined_checkm = "data/combined_checkm.tsv",
        combined_checkm_full = "data/combined_checkm_full.tsv",
    output:
        representative_checkm = "data/representative_checkm.tsv",
        representative_checkm_full = "data/representative_checkm_full.tsv",
    shell:
        "grep \'Bin Id\' {input.combined_checkm} >> {output.representative_checkm}; "
        "grep \'Bin Id\' {input.combined_checkm_full} >> {output.representative_checkm_full}; "
        "cut -f1 {input.dereplicated_clusters} | uniq | parallel -j 1 \"grep {{.}} {input.combined_checkm} >> {output.representative_checkm}\"; "
        "cut -f1 {input.dereplicated_clusters} | uniq | parallel -j 1 \"grep {{.}} {input.combined_checkm_full} >> {output.representative_checkm_full}\"; "


rule generate_pggb_input:
    input:
        dereplicated_clusters = "data/dereplicated_clusters.txt",
        # representative_info = "data/representative_info.txt"
    output:
        combined_mags_directory = directory("data/combined_mags/"),
        pggb_metadata = temp('data/pggb_metadata.tsv')
    run:
        from Bio import SeqIO
        import pandas as pd
        dereplicated_clusters = pd.read_csv(input.dereplicated_clusters, sep='\t', header=None)

        os.makedirs(output.combined_mags_directory, exist_ok=True)
        prev_cluster = ""
        taxonomy = ""
        sequences = []
        n_haplotypes = 0
        metadata_table = []
        for i in range(dereplicated_clusters.shape[0]):

            if prev_cluster != dereplicated_clusters.iloc[i, 0]:
                if len(sequences) != 0 and n_haplotypes > 1:
                    # print the fasta file
                    representative_sample_name = prev_cluster.split("/")[-4]
                    representative_name = prev_cluster.split("/")[-1].strip(".fna")
                    representative_name = f"{representative_sample_name}_{representative_name}"
                    output_filename = f"data/combined_mags/sequences_{representative_name}.fna"
                    SeqIO.write(sequences, output_filename, "fasta")
                    metadata_table.append([output_filename, n_haplotypes, representative_name])
                elif n_haplotypes == 1:
                    print(f"Cannot generate pangenome for single genome cluster: {dereplicated_clusters.iloc[i, 1]}")
                    print("Skipping...")
                # taxonomy = bin_info[bin_info['Bin Id'] == prev_cluster]['classification']
                prev_cluster = dereplicated_clusters.iloc[i, 0]
                sequences = []
                n_haplotypes = 0

            sample_name = prev_cluster.split("/")[-4]
            name = prev_cluster.split("/")[-1].strip(".fna")
            n_haplotypes += 1
            for record in SeqIO.parse(dereplicated_clusters.iloc[i, 1], "fasta"):
                record.id = f"{sample_name}_{name}_{record.id}"
                sequences.append(record)

        pd.DataFrame(metadata_table, columns=None).to_csv(output.pggb_metadata, sep='\t', index=False, header=None)

rule generate_pangenomes:
    input:
        pggb_metadata = 'data/pggb_metadata.tsv'
    output:
        pangenomes_output = directory("pangenomes/")
    params:
        derep_ani = config["ani"],
        pggb_params = config["pggb_params"]
    resources:
        mem_mb=int(config["max_memory"])*512
    conda:
        'envs/pggb.yaml'
    threads:
        config['max_threads']
    benchmark:
        'benchmarks/pggb.benchmark.txt'
    shell:
        'mkdir -p {output.pangenomes_output}; '
        'cat {input.pggb_metadata} | parallel -j1 --colsep \'\t\' '
        '"samtools faidx {{1}}; pggb -i {{1}} -n {{2}} -o {output.pangenomes_output}/{{3}} -t {threads} -p {params.derep_ani} -m {params.pggb_params}"'

rule complete_cluster:
    input:
        dereplicated_clusters = "data/dereplicated_clusters.txt",
        representative_checkm = "data/representative_checkm.tsv",
        representative_checkm_full = "data/representative_checkm_full.tsv",
        pggb_done = "benchmarks/pggb.benchmark.txt"
    output:
        representative_paths = "data/representative_paths.txt",
        representative_info = "data/representative_info.txt"
    params:
        previous_runs = config["previous_runs"]
    run:
        import pandas as pd

        clusters = pd.read_csv(input.dereplicated_clusters, sep="\t", header=None)
        representatives = clusters[0].apply(lambda bin: bin.replace(".fna", "")).unique()
        previous_runs = params.previous_runs
        # print(previous_runs)
        runs_info = []

        for run in previous_runs:
            print(f"Gathering: {run}")
            try:
                df_info = pd.read_csv(os.path.abspath(run) + "/bins/bin_info.tsv", sep="\t")
                df_info['Bin Id'] = df_info['Bin Id'].apply(lambda mag: os.path.abspath(run) + "/bins/final_bins/" + str(mag))
                runs_info.append(df_info)
            except (FileNotFoundError, IndexError):
                pass


        try:
            concat_bac = pd.concat(runs_info)
            concat_bac = concat_bac[concat_bac['Bin Id'].isin(representatives)]
            concat_bac.to_csv(output.representative_info, index=False, sep="\t")
        except ValueError:
            # No values to concat
            pass

        pd.DataFrame(representatives).to_csv("data/representative_paths.txt", index=False, header=False, sep="\t")

