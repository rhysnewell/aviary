import pandas as pd
from Bio import SeqIO
import os
from glob import glob

def find_circular(checkm_output, assembly_info_path):
    bin_column = "Name"

    circular_contigs = []
    circular_bps = []
    circular_fractions = []

    assembly_info = pd.read_csv(assembly_info_path, sep="\t")

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

def get_taxonomy(rename_columns="Name"):
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

    # checkm2 file for all bins
    checkm_output = pd.read_csv(snakemake.input.checkm2_done, sep='\t')
    coverage_file.rename({"Genome" : checkm_output.columns[0]}, inplace=True, axis=1)


    assembly_dir = f"data/{snakemake.config.get('long_read_assembler', 'myloasm')}"
    assembly_info_path = os.path.join(assembly_dir, "assembly_info.txt")

    if os.path.isfile(assembly_info_path):
        checkm_output = find_circular(checkm_output, assembly_info_path)

    taxa = get_taxonomy(checkm_output.columns[0])

    merged_out = pd.merge(checkm_output, coverage_file, on=[checkm_output.columns[0]], how="left")
    merged_out = pd.merge(merged_out, taxa, on=[merged_out.columns[0]], how="left")
    merged_out.to_csv(snakemake.output.bin_stats, sep='\t', index=False)

    checkm_minimal = checkm_output[[checkm_output.columns[0], "Completeness", "Contamination"]]

    checkm_minimal.to_csv(snakemake.output.checkm_minimal, sep="\t", index=False)
