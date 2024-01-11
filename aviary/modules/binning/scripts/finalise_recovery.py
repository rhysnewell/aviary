import os
from pathlib import Path
import glob

def check_and_remove_base_file(file_path):
    file_name = os.path.basename(coverage_file)
    if os.path.exists(file_name):
        os.remove(file_name)

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
        check_and_remove_base_file(coverage_file)
        os.symlink(f"../{coverage_file}", "./")
    
    if len(contig_coverage) > 0:
        check_and_remove_base_file(contig_coverage)
        os.symlink(f"../{contig_coverage}", "./")
    
    os.chdir('..')
    if len(gtdbtk) > 0:
        check_and_remove_base_file(output_taxonomy)
        os.symlink(f"{gtdbtk}", output_taxonomy)
    if len(singlem) > 0:
        check_and_remove_base_file(output_singlem)
        os.symlink(f"{singlem}", output_singlem)
    
    for f in glob.glob('data/binning_bams/*.ba*'):
        os.remove(f)

    Path('bins/done').touch()
