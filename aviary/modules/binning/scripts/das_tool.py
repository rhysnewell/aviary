import logging
import os
import sys

import extern

if __name__ == '__main__':
    unrefined_binners_to_use = [
        ('concoct', 'fa'),
        ('maxbin2', 'fasta'),
        ('vamb', 'fna'),
        ('comebin', 'fa'),
        ('taxvamb', 'fna'),
        ]
    refined_binners_to_use = [
        ('rosella', 'fna'),
        ('semibin', 'fna'),
        ]

    # N.B. specifying "metabat" will skip both MetaBAT1 and MetaBAT2.
    metabats = ['metabat_sspec', 'metabat_ssens', 'metabat_sens', 'metabat_spec']

    binners = []
    for (binner, extension) in unrefined_binners_to_use:
        if binner not in snakemake.config['skip_binners']:
            extra = ''
            if binner == 'vamb' or binner == 'taxvamb':
                extra = 'bins/'
            elif binner == 'comebin':
                extra = 'comebin_res/comebin_res_bins/'

            binners.append((f'{binner}_bins/'+extra, extension, f'data/{binner}_bins.tsv'))

    for (binner, extension) in refined_binners_to_use:
        if binner not in snakemake.config['skip_binners']:
            binners.append((binner+'_refined/final_bins/', extension, f'data/{binner}_refined_bins.tsv'))

    for metabat in metabats:
        if metabat not in snakemake.config['skip_binners']:
            binners.append((metabat.replace('metabat','metabat_bins'), 'fa', f'data/{metabat}_bins.tsv'))
    if 'metabat2' not in snakemake.config['skip_binners']:
        binners.append(('metabat2_refined/final_bins/', 'fna', f'data/metabat2_refined_bins.tsv'))

    logfile = snakemake.log[0]
    logging.basicConfig(filename=logfile, level=logging.INFO)
    logging.info("Using the following binners: " + str(binners))

    if len(binners) == 0:
        logging.error("All binners have been skipped, so DAS_tool cannot be run.")
        sys.exit(1)

    bin_definition_files = []
    for binner, extension, bin_definition_file in binners:
        extern.run(f'Fasta_to_Scaffolds2Bin.sh -i data/{binner} -e {extension} >{bin_definition_file}  2>> {logfile}')
        if os.path.getsize(bin_definition_file) == 0:
            logging.warning(f'Bin definition file {bin_definition_file} is empty, suggesting that {binner} failed or did not not create any output bins.')
        else:
            bin_definition_files.append(bin_definition_file)
    logging.info("Bin definition files created: " + str(bin_definition_files))

    scaffold2bin_files = ','.join(bin_definition_files)

    das_tool_command = f'DAS_Tool --search_engine diamond --write_bin_evals 1 --write_bins 1 -t {snakemake.threads} --score_threshold -42 \
        -i {scaffold2bin_files} \
        -c {snakemake.input.fasta} \
        -o data/das_tool_bins_pre_refine/das_tool >> {logfile} 2>&1'
    logging.info("Running DAS_Tool with command: " + das_tool_command)
    extern.run(das_tool_command)
