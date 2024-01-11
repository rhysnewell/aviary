import os
from pathlib import Path
import glob

if __name__ == '__main__':
    final_bins = snakemake.input.final_bins
    coverage_file = snakemake.input.coverm
    contig_coverage = snakemake.input.contig_coverage
    gtdbtk = snakemake.input.gtdbtk
    singlem = snakemake.input.singlem

    os.chdir('bins/')

    if len(coverage_file) > 0:
        os.symlink(f"../{coverage_file}", "./")
    
    if len(contig_coverage) > 0:
        os.symlink(f"../{contig_coverage}", "./")
    
    os.chdir('..')
    if len(gtdbtk) > 0:
        os.symlink(f"{gtdbtk}", "taxonomy")
    if len(singlem) > 0:
        os.symlink(f"{singlem}", "diversity")

    Path('bins/done').touch()

    for f in glob.glob('data/binning_bams/*.ba*'):
        os.remove(f)