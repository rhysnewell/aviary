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
        combined_checkm_file = temp("data/combined_checkm.tsv"),
        combined_checkm_stats = temp("data/combined_checkm_full.tsv"),
    params:
        previous_runs = list(config["previous_runs"])
    run:
        import pandas as pd

        previous_runs = params.previous_runs
        runs = [pd.read_csv(os.path.abspath(run) + "/bins/checkm_minimal.tsv", sep="\t") for run in previous_runs]
        infos = [pd.read_csv(os.path.abspath(run) + "/bins/bin_info.tsv", sep="\t") for run in previous_runs]

        for (df, run, info_df) in zip(runs, previous_runs, infos):
            df['Bin Id'] = df['Bin Id'].apply(lambda bin: os.path.abspath(run) + "/bins/final_bins/" + str(bin))
            info_df['Bin Id'] = info_df['Bin Id'].apply(lambda bin: os.path.abspath(run) + "/bins/final_bins/" + str(bin))

        concat = pd.concat(runs)
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
        "--min-completeness {params.min_completeness} --max-contamination {params.max_contamination} "
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

rule complete_cluster:
    input:
        dereplicated_clusters = "data/dereplicated_clusters.txt",
        representative_checkm = "data/representative_checkm.tsv",
        representative_checkm_full = "data/representative_checkm_full.tsv"
    output:
        representative_paths = "data/representative_paths.txt",
        # representative_taxa_bac = "data/representatives.bac120.summary.tsv",
        # representative_taxa_arc = "data/representatives.ar55.summary.tsv",
    params:
        previous_runs = config["previous_runs"]
    run:
        import pandas as pd
        import glob

        clusters = pd.read_csv(input.dereplicated_clusters, sep="\t", header=None)
        representatives = clusters[0].apply(lambda bin: bin.replace(".fna", "")).unique()
        previous_runs = params.previous_runs
        # print(previous_runs)
        runs_bac = []
        runs_arc = []

        for run in previous_runs:
            print(f"Gathering: {run}")
            try:
                df_bac = pd.read_csv(glob.glob(os.path.abspath(run) + "/taxonomy/gtdbtk.bac*.summary.tsv")[0], sep="\t")
                df_bac['user_genome'] = df_bac['user_genome'].apply(lambda mag: os.path.abspath(run) + "/bins/final_bins/" + str(mag))
                runs_bac.append(df_bac)
            except (FileNotFoundError, IndexError):
                pass

            try:
                df_arc = pd.read_csv(glob.glob(os.path.abspath(run) + "/taxonomy/gtdbtk.ar*.summary.tsv")[0], sep="\t")
                df_arc['user_genome'] = df_arc['user_genome'].apply(lambda mag: os.path.abspath(run) + "/bins/final_bins/" + str(mag))
                runs_arc.append(df_arc)
            except (FileNotFoundError, IndexError) as e:
                pass


        try:
            concat_bac = pd.concat(runs_bac)
            concat_bac = concat_bac[concat_bac['user_genome'].isin(representatives)]
            concat_bac.to_csv("data/representatives.bac120.summary.tsv", index=False, sep="\t")
        except ValueError:
            # No values to concat
            pass

        try:
            concat_arc = pd.concat(runs_arc)
            concat_arc = concat_arc[concat_arc['user_genome'].isin(representatives)]
            concat_arc.to_csv("data/representatives.ar55.summary.tsv", index=False, sep="\t")
        except ValueError:
            # No values to concate
            pass

        pd.DataFrame(representatives).to_csv("data/representative_paths.txt", index=False, header=False, sep="\t")

