import os
from pathlib import Path
import glob

def check_and_remove_base_file(file_path) -> str:
    file_name = os.path.basename(file_path)
    if os.path.exists(file_name):
        os.remove(file_name)

    return file_name


if __name__ == '__main__':
    final_bins = snakemake.input.final_bins
    coverage_file = snakemake.input.coverm
    contig_coverage = snakemake.input.contig_coverage
    gtdbtk = snakemake.input.gtdbtk
    singlem = snakemake.input.singlem

    output_taxonomy = "taxonomy"
    output_singlem = "diversity"

    os.chdir('bins/')

    if len(coverage_file) > 0:
        file_name = check_and_remove_base_file(coverage_file)
        os.symlink(f"../{coverage_file}", f"{file_name}")
    
    if len(contig_coverage) > 0:
        file_name = check_and_remove_base_file(contig_coverage)
        os.symlink(f"../{contig_coverage}", f"{file_name}")
    
    os.chdir('..')
    if len(gtdbtk) > 0:
        check_and_remove_base_file(output_taxonomy)
        os.symlink(f"{os.path.dirname(gtdbtk)}", output_taxonomy)
    if len(singlem) > 0:
        check_and_remove_base_file(output_singlem)
        os.symlink(f"{os.path.dirname(singlem)}", output_singlem)
    
    for f in glob.glob('data/binning_bams/*.ba*'):
        os.remove(f)

    Path('bins/done').touch()
