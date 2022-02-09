import pandas as pd
import subprocess
import shutil
import glob
import os


def refinery():
    """
    Main function that performs the refining process iteratively passing bins between
    rosella refine and checkm until either the max number of iterations is reached
    or there are no more contaminated bins
    """
    # These will change
    current_iteration = 0
    checkm_path = snakemake.input.checkm
    current_checkm = pd.read_csv(checkm_path, sep='\t', comment="[")
    bin_folder = "data/das_tool_bins/das_tool_DASTool_bins/"
    extension = "fa"

    # These will not change
    coverage = snakemake.input.coverage
    assembly = snakemake.input.fasta
    kmers = snakemake.input.kmers
    min_bin_size = snakemake.params.min_bin_size
    max_iterations = int(snakemake.params.max_iterations)
    pplacer_threads = int(snakemake.params.pplacer_threads)
    threads = int(snakemake.threads)
    max_contamination = int(snakemake.params.max_contamination)
    contaminated_bin_folder = "data/refined_bins/contaminated_bins"
    final_bins = "data/refined_bins/final_bins"

    os.makedirs(contaminated_bin_folder, exist_ok=True)
    os.makedirs(final_bins, exist_ok=True)

    final_checkm = current_checkm.copy().loc[current_checkm["Contamination"] <= 10].copy()
    final_checkm = move_finished_bins(final_checkm, bin_folder, extension, final_bins)


    while current_iteration < max_iterations:

        # delete previous contaminated set
        if current_iteration != 0:
            files = glob.glob(f'{contaminated_bin_folder}/*')
            [os.remove(file) for file in files]
            extension = "fna"

        # put contaminated bins in one folder. If no contaminated bins then break
        if not collect_contaminated_bins(
                current_checkm, max_contamination, bin_folder, extension, contaminated_bin_folder):
            current_checkm = pd.read_csv(f"{bin_folder}/checkm.out", sep='\t', comment='[')
            bins_to_keep = current_checkm[current_checkm["Contamination"] <= max_contamination]
            bins_to_keep = move_finished_bins(bins_to_keep, bin_folder, "fna", final_bins, current_iteration)
            final_checkm = pd.concat([final_checkm, bins_to_keep])
            break

        output_folder = f"data/refined_bins/rosella_refined_{current_iteration}"

        # Refine the contaminated bins
        refine(assembly, coverage, kmers, checkm_path,
               contaminated_bin_folder, extension, min_bin_size,
               threads, output_folder, max_contamination)

        # update the bin folder variable to the current refined folder
        bin_folder = output_folder

        # retrieve the checkm results for the refined bins
        get_checkm_results(bin_folder, threads, pplacer_threads)

        # update the checkm results and counter
        checkm_path = f"{bin_folder}/checkm.out"
        current_checkm = pd.read_csv(checkm_path, sep='\t', comment='[')

        bins_to_keep = current_checkm.copy().loc[current_checkm["Contamination"] <= max_contamination]
        bins_to_keep = move_finished_bins(bins_to_keep, bin_folder, "fna", final_bins, current_iteration)
        final_checkm = pd.concat([final_checkm, bins_to_keep])
        current_iteration += 1

    final_checkm.to_csv("data/checkm.out", sep='\t', index=False)
    final_output_folder = "bins/"
    os.makedirs(final_output_folder, exist_ok=True)
    if not os.path.exists(final_output_folder):
        os.symlink(os.path.abspath(final_bins), os.path.abspath(final_output_folder), target_is_directory=True)
    elif not os.path.islink(final_output_folder):
        os.rmdir(final_output_folder)
        os.symlink(os.path.abspath(final_bins), os.path.abspath(final_output_folder), target_is_directory=True)
    final_checkm.to_csv("bins/checkm.out", sep='\t', index=False)



def move_finished_bins(
        checkm,
        input_folder,
        input_extension,
        output_folder,
        current_iteration = None
):
    try:
        bins = checkm["Bin Id"]
        column_name = "Bin Id"
    except KeyError:
        # checkm2 input
        bins = checkm["Name"]
        column_name = "Name"

    new_bin_ids = []
    for bin_id in bins:

        if current_iteration is not None:
            new_bin_id = f"{bin_id}_{current_iteration}"
        else:
            new_bin_id = bin_id

        new_bin_ids.append(new_bin_id)
        try:
            shutil.copy(f"{input_folder}/{bin_id}.{input_extension}", f"{output_folder}/{new_bin_id}.fna")
        except FileNotFoundError:
            # File already moved
            continue

    checkm[column_name] = new_bin_ids

    return checkm


def collect_contaminated_bins(
        checkm,
        max_contamination,
        bin_folder,
        extension,
        output_folder
):
    contaminated = checkm[checkm["Contamination"] > max_contamination]

    # check that there is at least one contaminated bin
    if contaminated.shape[0] == 0: return False

    try:
        bins = contaminated["Bin Id"]
    except KeyError:
        # checkm2 input
        bins = contaminated["Name"]


    for bin_id in bins:
        try:
            shutil.copy(f"{bin_folder}/{bin_id}.{extension}", f"{output_folder}/{bin_id}.{extension}")
        except FileNotFoundError:
            # file already moved
            continue

    return True


def refine(
        assembly, coverage, kmers, checkm,
        bin_folder, extension, min_bin_size,
        threads, output_folder, max_contamination,
):

    subprocess.Popen(f"rosella refine -a {assembly} --coverage-values {coverage} --kmer-frequencies {kmers} "
                     f"-d {bin_folder} -x {extension} --checkm-file {checkm} --max-contamination {max_contamination} "
                     f"--min-bin-size {min_bin_size} -t {threads} -o {output_folder}", shell=True).wait()

def get_checkm_results(
        refined_folder,
        threads,
        pplacer_threads,
):
    subprocess.Popen(f"checkm lineage_wf -t {threads} --pplacer_threads {pplacer_threads} -x fna "
                     f"--tab_table -f {refined_folder}/checkm.out {refined_folder} {refined_folder}/checkm", shell=True).wait()

if __name__ == '__main__':
    refinery()