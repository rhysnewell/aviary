import pandas as pd
from Bio import SeqIO
import os
from glob import glob

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
        df_bac = pd.read_csv(glob("data/gtdbtk/gtdbtk.bac*.summary.tsv")[0], sep="\t")
        taxa.append(df_bac)
    except (FileNotFoundError, IndexError):
        pass

    try:
        df_arc = pd.read_csv(glob("data/gtdbtk/gtdbtk.ar*.summary.tsv")[0], sep="\t")
        taxa.append(df_arc)
    except (FileNotFoundError, IndexError) as e:
        pass
    
    try:
        taxa = pd.concat(taxa)
        taxa.rename({'user_genome' : rename_columns}, inplace=True, axis=1)
    except ValueError:
        taxa = pd.DataFrame(columns=[rename_columns])
    return taxa


if __name__ == "__main__":
    try:
        coverage_file = pd.read_csv(snakemake.input.coverage_file, sep='\t')
    except ValueError:
        coverage_file = pd.DataFrame(columns=["Genome"])

    # checkm file for all bins
    checkm1_output = pd.read_csv(snakemake.input.checkm1_done, sep='\t', comment="[")

    checkm2_output = pd.read_csv(snakemake.input.checkm2_done, sep='\t')

    checkm1_output.rename({'Completeness' : 'Completeness (CheckM1)', 'Contamination' : 'Contamination (CheckM1)'}, inplace=True, axis=1)
    checkm2_output.rename({'Name' : checkm1_output.columns[0], 'Completeness' : 'Completeness (CheckM2)', 'Contamination' : 'Contamination (CheckM2)'}, inplace=True, axis=1)

    checkm_output = pd.merge(checkm1_output, checkm2_output, on=[checkm1_output.columns[0]])
    is_checkm1 = "Bin Id" in checkm_output.columns
    coverage_file.rename({"Genome" : checkm_output.columns[0]}, inplace=True, axis=1)


    if os.path.isfile("data/flye/assembly_info.txt"):
        checkm_output = find_circular(checkm_output, is_checkm1)

    taxa = get_taxonomy(checkm_output.columns[0])

    merged_out = pd.merge(checkm_output, coverage_file, on=[checkm_output.columns[0]], how="left")
    merged_out = pd.merge(merged_out, taxa, on=[merged_out.columns[0]], how="left")
    merged_out.to_csv(snakemake.output.bin_stats, sep='\t', index=False)

    checkm_minimal = checkm_output[["Bin Id",  "Marker lineage",  "# genomes", "# markers", "# marker sets",
                                    "0", "1", "2", "3", "4", "5+", "Completeness (CheckM1)", "Contamination (CheckM1)",
                                    "Completeness (CheckM2)", "Contamination (CheckM2)", "Strain heterogeneity"]]

    checkm_minimal.to_csv(snakemake.output.checkm_minimal, sep="\t", index=False)