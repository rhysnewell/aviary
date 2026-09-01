#!/usr/bin/env python3
import logging
import os
import sys
import argparse
import subprocess  # replacing extern

def run(command):
    """Simple replacement for extern.run"""
    return subprocess.run(command, shell=True, check=True)

def reprefix_semibin_contig2bin(tsv_path: str, assembly_fasta: str) -> None:
    """Re-attach the 'sample:' prefix to SemiBin multi-mode bin contigs.

    In multi mode, SemiBin2 writes per-sample bin files (e.g.
    assembly2_SemiBin_0.fna) but STRIPS the 'sample:' prefix from the contig
    headers inside them, encoding the sample in the filename/bin-id instead.
    large_contigs.fasta keeps the prefix (assembly2:NODE_1), so the bare bin
    names (NODE_1) no longer match the assembly and DAS_Tool fails with
    "Contigs of contig2bin files not found in assembly". The sample survives in
    the bin id (<sample>_SemiBin_<n>), so re-derive it and restore the prefix.

    Self-guarding: only rewrites a contig when its bare name is ABSENT from the
    assembly and the prefixed candidate is PRESENT. Single-sample runs have bare
    names that already match, so they are left untouched.
    """
    assembly_names = {
        line[1:].split()[0]
        for line in open(assembly_fasta)
        if line.startswith('>')
    }
    changed = 0
    fixed_lines = []
    for line in open(tsv_path):
        if not line.strip():
            continue
        contig, binid = line.rstrip('\n').split('\t')
        if contig not in assembly_names and '_SemiBin_' in binid:
            candidate = f"{binid.rsplit('_SemiBin_', 1)[0]}:{contig}"
            if candidate in assembly_names:
                contig = candidate
                changed += 1
        fixed_lines.append(f'{contig}\t{binid}')
    with open(tsv_path, 'w') as f:
        f.write('\n'.join(fixed_lines) + '\n')
    logging.info(f'semibin contig2bin: re-prefixed {changed} contigs to match assembly')

def parse_args():
    parser = argparse.ArgumentParser(description='Wrapper script for DAS Tool to process multiple binners')
    
    # Input files
    parser.add_argument('--skip-binners', nargs='*', required=True, help='List of binners to skip')
    parser.add_argument('--fasta', required=True, help='Path to the assembly fasta file')
    
    # Process parameters
    parser.add_argument('--threads', required=True, type=int, help='Number of threads to use')
    
    # Log parameter
    parser.add_argument('--log', required=True, help='Path to log file')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    
    unrefined_binners_to_use = [
        ('concoct', 'fa'),
        ('maxbin2', 'fasta'),
        ('vamb', 'fna'),
        ('comebin', 'fa'),
        ('taxvamb', 'fna'),
        ('quickbin', 'fna'),
    ]
    refined_binners_to_use = [
        ('rosella', 'fna'),
        ('semibin', 'fna'),
    ]

    # N.B. specifying "metabat" will skip both MetaBAT1 and MetaBAT2.
    metabats = ['metabat_sspec', 'metabat_ssens', 'metabat_sens', 'metabat_spec']

    binners = []
    for (binner, extension) in unrefined_binners_to_use:
        if binner not in args.skip_binners:
            extra = ''
            if binner == 'vamb' or binner == 'taxvamb':
                extra = 'bins/'
            elif binner == 'comebin':
                extra = 'comebin_res/comebin_res_bins/'

            binners.append((f'{binner}_bins/'+extra, extension, f'data/{binner}_bins.tsv'))

    for (binner, extension) in refined_binners_to_use:
        if binner not in args.skip_binners:
            binners.append((binner+'_refined/final_bins/', extension, f'data/{binner}_refined_bins.tsv'))

    for metabat in metabats:
        if metabat not in args.skip_binners:
            binners.append((metabat.replace('metabat','metabat_bins'), 'fa', f'data/{metabat}_bins.tsv'))
    if 'metabat2' not in args.skip_binners:
        binners.append(('metabat2_refined/final_bins/', 'fna', f'data/metabat2_refined_bins.tsv'))

    logfile = args.log
    logging.basicConfig(filename=logfile, level=logging.INFO)
    logging.info("Using the following binners: " + str(binners))

    if len(binners) == 0:
        logging.error("All binners have been skipped, so DAS_tool cannot be run.")
        sys.exit(1)

    bin_definition_files = []
    for binner, extension, bin_definition_file in binners:
        # metabat2 adds \t total_depth=X \t sample_depths=X to contig names, so we need to remove those.
        # Causes dastool error: ERROR: Cannot read scaffold2bin file: data/metabat2_refined_bins.tsv.	Please check file/format. Should be: scaffold_name \t bin_id 
        if 'metabat2' in binner:
            fix_metabat_cmd = "| awk -F'\t' '{print $1 \"\t\" $NF}'"
        else:
            fix_metabat_cmd = ""

        run(f'Fasta_to_Contig2Bin.sh -i data/{binner} -e {extension} {fix_metabat_cmd} >{bin_definition_file}  2>> {logfile}')
        if binner.startswith('semibin') and os.path.getsize(bin_definition_file) > 0:
            reprefix_semibin_contig2bin(bin_definition_file, args.fasta)
        if os.path.getsize(bin_definition_file) == 0:
            logging.warning(f'Bin definition file {bin_definition_file} is empty, suggesting that {binner} failed or did not not create any output bins.')
        else:
            bin_definition_files.append(bin_definition_file)

    if len(bin_definition_files) == 0:
        logging.error("No bins were found, so DAS_tool cannot be run.")
        raise RuntimeError("No bins were found, so DAS_tool cannot be run.")

    logging.info("Bin definition files created: " + str(bin_definition_files))


    scaffold2bin_files = ','.join(bin_definition_files)

    das_tool_command = f'DAS_Tool --search_engine diamond --write_bin_evals --write_bins -t {args.threads} --score_threshold -42 \
        -i {scaffold2bin_files} \
        -c {args.fasta} \
        -o data/das_tool_bins_pre_refine/das_tool >> {logfile} 2>&1'
    logging.info("Running DAS_Tool with command: " + das_tool_command)
    run(das_tool_command)
