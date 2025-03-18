import pandas as pd
from subprocess import run, STDOUT
import shutil
import glob
import os
import pandas.errors as e
import datetime

def refinery():
    """
    Main function that performs the refining process iteratively passing bins between
    rosella refine and checkm until either the max number of iterations is reached
    or there are no more contaminated bins
    """
    # These will change
    current_iteration = 0
    checkm_path = snakemake.input.checkm

    # These will not change
    coverage = snakemake.input.coverage
    assembly = snakemake.input.fasta
    final_refining = snakemake.params.final_refining
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    # kmers = snakemake.input.kmers
    if os.path.isfile("data/rosella_bins/kmer_frequencies.tsv"):
        kmers = "data/rosella_bins/kmer_frequencies.tsv"
    else:
        kmers = None
    min_bin_size = snakemake.params.min_bin_size
    max_iterations = int(snakemake.params.max_iterations)
    max_retries = int(snakemake.params.max_retries)
    pplacer_threads = int(snakemake.params.pplacer_threads)
    threads = int(snakemake.threads)
    max_contamination = int(snakemake.params.max_contamination)
    contaminated_bin_folder = snakemake.params.output_folder + "/contaminated_bins"
    final_bins = snakemake.params.output_folder + "/final_bins"
    bin_folder = snakemake.params.bin_folder
    extension = snakemake.params.extension
    bin_prefix = snakemake.params.bin_prefix


    try:
        os.makedirs(contaminated_bin_folder)
    except FileExistsError:
        shutil.rmtree(contaminated_bin_folder)
        os.makedirs(contaminated_bin_folder)

    try:
        os.makedirs(final_bins)
    except FileExistsError:
        shutil.rmtree(final_bins)
        os.makedirs(final_bins)

    try:
        current_checkm = pd.read_csv(checkm_path, sep='\t', comment="[")
    except e.EmptyDataError:
        with open(log, "a") as logf:
            logf.write("No bins found in checkm file\n")
            logf.write("Skipping refinement\n")
        if max_iterations == 0:
            for bin in os.listdir(bin_folder):
                if bin.endswith(extension):
                    shutil.copy(f"{bin_folder}/{bin}", f"{final_bins}/{os.path.splitext(bin)[0]}.fna")

        open(f"{snakemake.params.output_folder}/done", "a").close()

        return

    final_checkm = current_checkm.copy().loc[current_checkm["Contamination"] <= snakemake.params.max_contamination].copy()
    final_checkm = move_finished_bins(final_checkm, bin_folder, extension, final_bins)
    previous_checkm_path = None
    unchanged_bins = None

    while current_iteration < max_iterations:
        with open(log, "a") as logf:
            # write white space for legibility
            logf.write("\n")
            logf.write("\n")

            logf.write(f"INFO: {datetime.datetime.now().strftime('%H:%M:%S')} - Refining iteration {current_iteration}\n")

        # delete previous contaminated set
        if current_iteration != 0:
            files = glob.glob(f'{contaminated_bin_folder}/*')
            [os.remove(file) for file in files]
            extension = "fna"

        # put contaminated bins in one folder. If no contaminated bins then break
        if not collect_contaminated_bins(
                current_checkm, max_contamination, bin_folder, extension, contaminated_bin_folder):
            break

        output_folder = f"{snakemake.params.output_folder}/rosella_refined_{current_iteration + 1}"

        if os.path.exists(output_folder):
            shutil.rmtree(output_folder)

        with open(log, "a") as logf:
            # write white space for legibility
            logf.write("\n")
            logf.write(f"INFO: {datetime.datetime.now().strftime('%H:%M:%S')} - Rosella refine iteration {current_iteration}\n")
        
        # Refine the contaminated bins
        kmers = refine(assembly, coverage, kmers, checkm_path,
                   contaminated_bin_folder, extension, min_bin_size,
                   threads, output_folder, max_contamination, max_retries, f"{bin_prefix}_refined_{current_iteration + 1}", log)

        # update the bin folder variable to the current refined folder
        bin_folder = f"{output_folder}/refined_bins/"
        unchanged_bins = f"{output_folder}/unchanged_bins/"

        # move the unchanged bins to the final bins folder
        if os.path.exists(unchanged_bins):
            for mag in os.listdir(unchanged_bins):
                shutil.move(f"{unchanged_bins}/{mag}", f"{final_bins}/{os.path.splitext(mag)[0]}.fna")

        # retrieve the checkm results for the refined bins
        try:
            with open(log, "a") as logf:
                # write white space for legibility
                logf.write("\n")
                logf.write(f"INFO: {datetime.datetime.now().strftime('%H:%M:%S')} - CheckM iteration {current_iteration}\n")
            # count how many bins in bin_folder, bins end in 'fna'
            bin_count = 0
            for bin_file in os.listdir(bin_folder):
                if bin_file.endswith("fna"):
                    bin_count += 1
            
            if bin_count == 0:
                with open(log, "a") as logf:
                    logf.write("No bins to refine\n")
                    logf.write("Skipping refinement\n")
                # No bins to refine, break out and move on
                break
            else:
                with open(log, "a") as logf:
                    logf.write(f"Running CheckM on {bin_count} bins\n")

            get_checkm_results(bin_folder, threads, pplacer_threads, log, final_refining)
            # update the checkm results and counter
            checkm_path = f"{bin_folder}/checkm.out"
            current_checkm = pd.read_csv(checkm_path, sep='\t', comment='[')

            bins_to_keep = current_checkm.copy().loc[current_checkm["Contamination"] <= max_contamination]
            bins_to_keep = move_finished_bins(bins_to_keep, bin_folder, "fna", final_bins, None) # None because bin tag handles refined bin name
            final_checkm = pd.concat([final_checkm, bins_to_keep])

            # also handle unchanged bins
            if os.path.exists(unchanged_bins):
                # list the bins in unchanged_bins
                if previous_checkm_path is None:
                    previous_checkm_path = snakemake.input.checkm
                
                bin_id_accessor_str = "Bin Id"

                previous_checkm = pd.read_csv(previous_checkm_path, sep='\t', comment='[')
                if "Bin Id" not in previous_checkm.columns:
                    # checkm2 input
                    bin_id_accessor_str = "Name"

                bin_ids = []
                for bin_path in os.listdir(unchanged_bins):
                    # copy the bin to the final bins folder
                    bin_id = os.path.splitext(bin_path)[0]
                    bin_ids.append(bin_id)
                    shutil.copy(f"{unchanged_bins}/{bin}", f"{final_bins}/{bin_id}.fna")
                    # add the bin to the final checkm
                    # row = previous_checkm.copy().loc[previous_checkm["Bin Id"] == bin_id]
                    # final_checkm = pd.concat([final_checkm, row])
                final_checkm = pd.concat([final_checkm, previous_checkm.copy().loc[previous_checkm[bin_id_accessor_str].isin(bin_ids)]])

            previous_checkm_path = checkm_path

            with open(log, "a") as logf:
                logf.write(f"INFO: {datetime.datetime.now().strftime('%H:%M:%S')} - Ending iteration: {current_iteration}\n")

            current_iteration += 1
        except FileNotFoundError:
            with open(log, "a") as logf:
                logf.write("No bins to refine\n")
                logf.write(f"Ending refine on iteration: {current_iteration}\n")
            # No bins to refine, break out and move on
            break
    

    if final_refining:
        final_checkm.to_csv("data/checkm.out", sep='\t', index=False)
        final_output_folder = "bins/final_bins"
        if os.path.exists("bins/"): # remove pre-existing output, preventing bin duplication
            shutil.rmtree("bins/")
        os.makedirs("bins/")
        if not os.path.exists(final_output_folder):
            os.symlink("../" + final_bins, final_output_folder, target_is_directory=True)
        elif not os.path.islink(final_output_folder):
            os.rmdir(final_output_folder)
            os.symlink("../" + final_bins, final_output_folder, target_is_directory=True)
        final_checkm.to_csv("bins/checkm.out", sep='\t', index=False)
    else:
        with open(log, "a") as logf:
            logf.write("Refinery finished.\n")
        
        # final_checkm.to_csv(f"{snakemake.params.output_folder}/intermediate_checkm.out", sep='\t', index=False)
        open(f"{snakemake.params.output_folder}/done", "a").close()
    
    # in output folder, we then delete directories starting with rosella_refined_
    for folder in os.listdir(snakemake.params.output_folder):
        if folder.startswith("rosella_refined_"):
            shutil.rmtree(f"{snakemake.params.output_folder}/{folder}")



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
        # print(f"{input_folder}/{bin_id}")
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
    assembly,
    coverage,
    kmers,
    checkm,
    bin_folder,
    extension,
    min_bin_size,
    threads,
    output_folder,
    max_contamination,
    max_retries,
    bin_tag,
    log,
):
    if kmers is None:
        rosella_cmd = f"rosella refine -r {assembly} -C {coverage} -d {bin_folder} -x {extension} --checkm-results {checkm} --max-contamination {max_contamination} --min-bin-size {min_bin_size} -t {threads} -o {output_folder} --bin-tag {bin_tag}".split()
        with open(log, "a") as logf:
            run(rosella_cmd, stdout=logf, stderr=STDOUT)

        os.makedirs("data/rosella_bins/", exist_ok=True)
        shutil.copyfile(f"{output_folder}/kmer_frequencies.tsv", f"data/rosella_bins/kmer_frequencies.tsv")
        kmers = "data/rosella_bins/kmer_frequencies.tsv"
    else:
        rosella_cmd = f"rosella refine -r {assembly} -C {coverage} --kmer-frequency-file {kmers} -d {bin_folder} -x {extension} --checkm-results {checkm} --max-contamination {max_contamination} --min-bin-size {min_bin_size} -t {threads} -o {output_folder} --bin-tag {bin_tag}".split()

        with open(log, "a") as logf:
            run(rosella_cmd, stdout=logf, stderr=STDOUT)

    return kmers

def get_checkm_results(
        refined_folder,
        threads,
        pplacer_threads,
        log,
        final_refining=False,
):

    checkm_cmd = f"checkm lineage_wf -t {threads} --pplacer_threads {pplacer_threads} -x fna --tab_table -f {refined_folder}/checkm.out {refined_folder} {refined_folder}/checkm".split()
    with open(log, "a") as logf:
        run(checkm_cmd, stdout=logf, stderr=STDOUT)

    if final_refining:
        checkm_qa_cmd = f"checkm qa -o 2 --tab_table -f {refined_folder}/checkm.out {refined_folder}/checkm/lineage.ms {refined_folder}/checkm/".split()
        with open(log, "a") as logf:
            run(checkm_qa_cmd, stdout=logf, stderr=STDOUT)


if __name__ == '__main__':
    refinery()