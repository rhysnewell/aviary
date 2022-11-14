#!/usr/bin/env python
###############################################################################
# aviary.py - Info about aviary.py
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
import aviary.config.config as Config
from aviary.modules.processor import Processor, process_batch
from .__init__ import __version__
__author__ = "Rhys Newell"
__copyright__ = "Copyright 2022"
__credits__ = ["Rhys Newell"]
__license__ = "GPL3"
__maintainer__ = "Rhys Newell"
__email__ = "rhys.newell near hdr.qut.edu.au"
__status__ = "Development"

###############################################################################
# System imports
import sys
import argparse
import logging
import os
from datetime import datetime
import tempfile

# Debug
debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################
############################### - Exceptions - ################################

class BadTreeFileException(Exception):
    pass

###############################################################################
################################ - Functions - ################################
def centerify(text, width=-1):
  lines = text.split('\n')
  width = max(map(len, lines)) if width == -1 else width
  return '\n'.join(line.center(width) for line in lines)


def phelp():
    print(
"""

                    ......:::::: AVIARY ::::::......

           A comprehensive metagenomics bioinformatics pipeline

Metagenome assembly, binning, and annotation:
        assemble  - Perform hybrid assembly using short and long reads, 
                    or assembly using only short reads
        recover   - Recover MAGs from provided assembly using a variety 
                    of binning algorithms 
        annotate  - Annotate MAGs using EggNOG and GTBD-tk
        diversity - Perform strain diversity analysis of MAGs using Lorikeet
        complete  - Runs each stage of the pipeline: assemble, recover, 
                    annotate, diversity in that order.
        cluster   - Combines and dereplicates the MAGs from multiple Aviary runs
                    using Galah
        batch     - Run Aviary using a given workflow on a supplied batch of samples
                    and cluster the end result.

Isolate assembly, binning, and annotation:
        isolate   - Perform isolate assembly **PARTIALLY COMPLETED**
        
Utility modules:
        configure - Set or overwrite the environment variables for future runs.

"""
)


def str2bool(v):
    if isinstance(v, bool):
        return(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return(True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return(False)
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        phelp()
        return

    # Source the conda environment variables in case users have previously set
    # the variables using config but have not restarted the environment.
    try:
        Config.source_conda_env()
    except FileNotFoundError:
        Config.source_bashrc()

    ############################ ~ Main Parser ~ ##############################
    main_parser = argparse.ArgumentParser(prog='aviary',
                                          formatter_class=CustomHelpFormatter,
                                          add_help=False)
    main_parser.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')
    main_parser.add_argument('--verbosity',
                             help='1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging)',
                             type=int,
                             default=4)
    main_parser.add_argument('--log',
                             help='Output logging information to file',
                             default=False)
    subparsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    #~#~#~#~#~#~#~#~#~#~#~#~#~ Command groups ~#~#~#~#~#~#~#~#~#~#~#~#~#

    ####################################################################

    base_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                         add_help=False)

    base_group.add_argument(
        '-t', '--max-threads', '--max_threads',
        help='Maximum number of threads given to any particular process',
        dest='max_threads',
        default=8,
    )

    base_group.add_argument(
        '-p', '--pplacer-threads', '--pplacer_threads',
        help=argparse.SUPPRESS,
        dest='pplacer_threads',
        default=8,
    )

    base_group.add_argument(
        '-n', '--n-cores', '--n_cores',
        help='Maximum number of cores available for use. Must be >= to max_threads',
        dest='n_cores',
        default=16,
    )

    base_group.add_argument(
        '-m', '--max-memory', '--max_memory',
        help='Maximum memory for available usage in Gigabytes',
        dest='max_memory',
        default=250,
    )

    base_group.add_argument(
        '-o', '--output',
        help='Output directory',
        dest='output',
        default='./',
    )

    base_group.add_argument(
        '--conda-prefix', '--conda_prefix',
        help='Path to the location of installed conda environments, or where to install new environments. \n'
             'Can be configured within the `configure` subcommand',
        dest='conda_prefix',
        default=None,
    )

    base_group.add_argument(
        '--tmpdir', '--tempdir', '--tmp-dir', '--tmp_dir', '--tmp', '--temp', '--temp-dir', '--temp_dir',
        help='Path to the location that will be treated used for temporary files. If none is specified, the TMPDIR \n'
             'environment variable will be used. Can be configured within the `configure` subcommand',
        dest='tmpdir',
        default=tempfile.gettempdir(),
    )

    base_group.add_argument(
        '--default-resources',
        help='Snakemake resources used as is found at: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=resources#standard-resources \n'
             'NOTE: tmpdir is handled by the `tmpdir` command line parameter. ',
        dest='resources',
        default=" "
    )

    base_group.add_argument(
        '--dry-run', '--dry_run', '--dryrun',
        help='Perform snakemake dry run, tests workflow order and conda environments',
        type=str2bool,
        nargs='?',
        const=True,
        dest='dryrun',
        default=False,
    )

    base_group.add_argument(
        '--conda-frontend', '--conda_frontend',
        help='Which conda frontend to use, mamba is faster but harder to debug. Switch this to conda \n'
             'If experiencing problems installing environments',
        dest='conda_frontend',
        default="mamba",
        choices=["conda", "mamba"],
    )

    base_group.add_argument(
        '--clean',
        help='Clean up all temporary files. This will remove most BAM files and any FASTQ files \n'
             'generated from read filtering. Setting this to False is the equivalent of the --notemp \n'
             'option in snakemake. Useful for when running only part of a workflow as it avoids \n'
             'deleting files that would likely be needed in later parts of the workflow. \n'
             'NOTE: Not cleaning makes reruns faster but will incur the wrath of your sysadmin',
        type=str2bool,
        nargs='?',
        const=True,
        dest='clean',
        default=True,
    )

    base_group.add_argument(
        '--build',
        help='Build conda environments and then exits. Equivalent to \"--snakemake-cmds \'--conda-create-envs-only True \' \"',
        type=str2bool,
        nargs='?',
        const=True,
        dest='build',
    )

    base_group.add_argument(
        '--download', '--download',
        help='Downloads the required GTDB, EggNOG, & CheckM2 databases if required',
        type=str2bool,
        nargs='?',
        const=True,
        dest='download',
    )

    base_group.add_argument(
        '--rerun-triggers', '--rerun_triggers',
        help='Specify which kinds of modifications will trigger rules to rerun',\
        dest='rerun_triggers',
        default="mtime",
        nargs="*",
        choices=["mtime","params","input","software-env","code"]
    )

    base_group.add_argument(
        '--snakemake-cmds',
        help='Additional commands to supplied to snakemake in the form of a single string '
             'e.g. "--print-compilation True". \n '
             'NOTE: Most commands in snakemake -h are valid but some commands may clash with commands \n '
             'aviary directly supplies to snakemake. Please make sure your additional commands don\'t clash.',
        dest='cmds',
        default='',
    )

    ####################################################################
    qc_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                               add_help=False)

    qc_group.add_argument(
        '-g', '--gold-standard-assembly', '--gold_standard_assembly',
        help='Gold standard assembly to compare either the Aviary assembly or a given input assembly against',
        dest="gold_standard",
        default=['none'],
        nargs='*'
    )

    qc_group.add_argument(
        '--gsa-mappings', '--gsa_mappings',
        help='CAMI I & II GSA mappings',
        dest="gsa_mappings",
        default='none'
    )

    qc_group.add_argument(
        '-r', '--reference-filter', '--reference_filter',
        help='Reference filter file to aid in the assembly',
        dest="reference_filter",
        nargs=1,
        default='none'
    )

    qc_group.add_argument(
        '--min-read-size', '--min_read_size',
        help='Minimum long read size when filtering using Filtlong',
        dest="min_read_size",
        default=250
    )

    qc_group.add_argument(
        '--min-mean-q', '--min_mean_q',
        help='Minimum mean quality threshold',
        dest="min_mean_q",
        default=50
    )

    qc_group.add_argument(
        '--keep-percent', '--keep_percent',
        help='Percentage of reads passing quality thresholds kept by filtlong',
        dest="keep_percent",
        default=100
    )


    ####################################################################

    short_read_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                         add_help=False)
    read_group_exclusive = short_read_group.add_mutually_exclusive_group()

    read_group_exclusive.add_argument(
        '-1', '--pe-1', '--paired-reads-1', '--paired_reads_1', '--pe1',
        help='A space separated list of forwards read files \n'
             'NOTE: If performing assembly and multiple files and longreads \n'
             '      are provided then only the first file will be used for assembly. \n'
             '      If no longreads are provided then all samples will be co-assembled \n'
             '      with megahit or metaspades depending on the --coassemble parameter\n',
        dest='pe1',
        nargs='*',
        default="none"
    )

    short_read_group.add_argument(
        '-2', '--pe-2', '--paired-reads-2', '--paired_reads_2', '--pe2',
        help='A space separated list of reverse read files \n'
             'NOTE: If performing assembly and multiple files and longreads \n'
             '      are provided then only the first file will be used for assembly. \n'
             '      If no longreads are provided then all samples will be co-assembled \n'
             '      with megahit or metaspades depending on the --coassemble parameter',
        dest='pe2',
        nargs='*',
        default="none"
    )

    read_group_exclusive.add_argument(
        '-i','--interleaved',
        help='A space separated list of interleaved read files \n'
             'NOTE: If performing assembly and multiple files and longreads \n'
             '      are provided then only the first file will be used for assembly. \n'
             '      If no longreads are provided then all samples will be co-assembled \n'
             '      with megahit or metaspades depending on the --coassemble parameter',
        dest='interleaved',
        nargs='*',
        default="none"
    )

    read_group_exclusive.add_argument(
        '-c', '--coupled',
        help='Forward and reverse read files in a coupled space separated list. \n'
             'NOTE: If performing assembly and multiple files and longreads \n'
             '      are provided then only the first file will be used for assembly. \n'
             '      If no longreads are provided then all samples will be co-assembled \n'
             '      with megahit or metaspades depending on the --coassemble parameter',
        dest='coupled',
        nargs='*',
        default="none"
    )

    read_group_exclusive.add_argument(
        '--min-percent-read-identity-short', '--min_percent_read_identity_short',
        help='Minimum percent read identity used by CoverM for short-reads \n'
             'when calculating genome abundances.',
        dest='short_percent_identity',
        default='95'
    )

    ####################################################################

    long_read_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                              add_help=False)
    long_read_group.add_argument(
        '-l', '--longreads', '--long-reads', '--long_reads',
        help='A space separated list of long-read read files. '
             'NOTE: If performing assembly and multiple long read files are provided, \n'
             '      then only the first file is used for assembly. This behaviour might change in future.',
        dest='longreads',
        nargs='*',
        default="none"
    )

    long_read_group.add_argument(
        '-z', '--longread-type', '--longread_type', '--long_read_type', '--long-read-type',
        help='Whether the sequencing platform and technology for the longreads. \n'
             '"rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS \n'
             'reads, "ont" for Oxford Nanopore and "ont_hq" for Oxford Nanopore high quality reads (Guppy5+ or Q20) \n',
        dest='longread_type',
        nargs=1,
        default="ont",
        choices=["ont","ont_hq", "rs", "sq", "ccs"],
    )

    long_read_group.add_argument(
        '--min-percent-read-identity-long', '--min_percent_read_identity_long',
        help='Minimum percent read identity used by CoverM for long-reads'
             'when calculating genome abundances.',
        dest='long_percent_identity',
        default='85'
    )

    ####################################################################

    annotation_group = argparse.ArgumentParser(add_help=False)

    # annotation_group.add_argument(
    #     '--enrichm-db-path', '--enrichm_db_path',
    #     help='Path to the local EnrichM Database files',
    #     dest='enrichm_db_path',
    #     default=Config.get_software_db_path('ENRICHM_DB', '--enrichm-db-path'),
    # )

    annotation_group.add_argument(
        '--gtdb-path', '--gtdb_path',
        help='Path to the local gtdb database files',
        dest='gtdb_path',
        default=None,
    )

    annotation_group.add_argument(
        '--eggnog-db-path', '--eggnog_db_path',
        help='Path to the local eggnog database files',
        dest='eggnog_db_path',
        default=None,
    )

    annotation_group.add_argument(
        '--checkm2-db-path', '--checkm2_db_path',
        help='Path to Checkm2 Database',
        dest='checkm2_db_path',
        required=False,
    )

    ####################################################################

    binning_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                            add_help=False)

    binning_group.add_argument(
        '-s', '--min-contig-size', '--min_contig_size',
        help='Minimum contig size in base pairs to be considered for binning',
        dest='min_contig_size',
        default=1500
    )

    binning_group.add_argument(
        '-b', '--min-bin-size', '--min_bin_size',
        help='Minimum bin size in base pairs for a MAG',
        dest='min_bin_size',
        default=200000
    )

    binning_group.add_argument(
        '--semibin-model', '--semibin_model',
        help='The environment model to passed to SemiBin. Can be one of: \n'
             'human_gut, dog_gut, ocean, soil, cat_gut, human_oral, mouse_gut, pig_gut, built_environment, wastewater, global',
        dest='semibin_model',
        default='global'
    )

    binning_group.add_argument(
        '--skip-binners', '--skip_binners', '--skip_binner', '--skip-binner',
        help='Optional list of binning algorithms to skip. Can be any combination of: \n'
             'rosella, semibin, metabat1, metabat2, metabat, vamb, concoct, maxbin2, maxbin \n'
             'Capitals will be auto-corrected. N.B. specifying "metabat" will skip both \n'
             'MetaBAT1 and MetaBAT2.',
        dest='skip_binners',
        nargs='*'
        # default=["maxbin2"]
    )

    ####################################################################
    mag_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                        add_help=False)
    mag_group_exclusive = mag_group.add_mutually_exclusive_group()

    # mag_group_exclusive.add_argument(
    #     '-f', '--genome-fasta-files', '--genome_fasta_files',
    #     help='MAGs to be annotated',
    #     dest='mags',
    #     nargs='*',
    #     required=False,
    # )

    mag_group_exclusive.add_argument(
        '-d', '--genome-fasta-directory', '--genome_fasta_directory',
        help='Directory containing MAGs to be annotated',
        dest='directory',
        required=False,
    )

    mag_group.add_argument(
        '-x', '--fasta-extension', '--fasta_extension',
        help='File extension of fasta files in --genome-fasta-directory',
        dest='ext',
        required=False,
        default='fna'
    )

    #####################################################################
    isolate_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                            add_help=False)

    isolate_group.add_argument(
        '--guppy-model', '--guppy_model',
        help='The guppy model used by medaka to perform polishing',
        dest='guppy_model',
        nargs=1,
        required=False,
        default='r941_min_high_g360'
    )

    isolate_group.add_argument(
        '--genome-size', '--genome_size',
        help='Approximate size of the isolate genome to be assembled',
        dest='genome_size',
        nargs=1,
        required=False,
        default=5000000
    )

    #####################################################################
    cluster_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter, add_help=False)

    cluster_group.add_argument(
        '--ani',
        help='Overall ANI level to dereplicate at with FastANI.',
        dest='ani',
        default=97
    )

    cluster_group.add_argument(
        '--precluster-ani', '--precluster_ani',
        help='Require at least this dashing-derived ANI for preclustering and to avoid FastANI on distant lineages within preclusters.',
        dest='precluster_ani',
        default=95
    )

    cluster_group.add_argument(
        '--precluster-method', '--precluster_method',
        help="method of calculating rough ANI for dereplication. 'dashing' for HyperLogLog, 'finch' for finch MinHash.",
        dest='precluster_method',
        default='dashing'
    )

    cluster_group.add_argument(
        '--min-completeness', '--min_completeness',
        help="Ignore genomes with less completeness than this percentage.",
        dest='min_completeness',
        default='none'
    )

    cluster_group.add_argument(
        '--max-contamination', '--max_contamination',
        help="Ignore genomes with more contamination than this percentage.",
        dest='max_contamination',
        default='none'
    )

    cluster_group.add_argument(
        '--use-checkm2-scores', '--use_checkm2_scores',
        help="Use CheckM2 completeness and contamination scores (if available) to perform Galah dereplication",
        type=str2bool,
        nargs='?',
        const=True,
        dest='use_checkm2_scores',
        default=False
    )

    cluster_group.add_argument(
        '--pggb-params', '--pggb_params',
        help="Parameters to be used with pggb, must be surrounded by quotation marks e.g. \'\'",
        dest='pggb_params',
        default='-k 79 -G 7919,8069'
    )

    #####################################################################
    # viral_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
    #                                         add_help=False)
    #
    # viral_group.add_argument(
    #     '--virsorter-data', '--virsorter_data',
    #     help='The guppy model used by medaka to perform polishing',
    #     dest='guppy_model',
    #     nargs=1,
    #     required=False,
    #     default='r941_min_high_g360'
    # )
    #
    # viral_group.add_argument(
    #     '--genome-size', '--genome_size',
    #     help='Approximate size of the isolate genome to be assembled',
    #     dest='genome_size',
    #     nargs=1,
    #     required=False,
    #     default=5000000
    # )

    assemble_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter, add_help=False)
    assemble_group.add_argument(
        '--use-unicycler', '--use_unicycler',
        help='Use Unicycler to re-assemble the metaSPAdes hybrid assembly. Not recommended for complex metagenomes.',
        type=str2bool,
        nargs='?',
        const=True,
        dest='use_unicycler',
        default=False,
    )

    assemble_group.add_argument(
        '--use-megahit', '--use_megahit',
        help='Specifies whether or not to use megahit if multiple for short-read only assembly',
        type=str2bool,
        nargs='?',
        const=True,
        dest='use_megahit',
        default=False,
    )

    assemble_group.add_argument(
        '--coassemble', '--co-assemble', '--co_assemble',
        help='Specifies whether or not, when given multiple input reads, to coassemble them. \n'
             'If False, Aviary will use the first set of short reads and first set of long reads to perform assembly \n'
             'All read files will still be used during the MAG recovery process for differential coverage.',
        type=str2bool,
        nargs='?',
        const=True,
        dest='coassemble',
        default=True,
    )

    assemble_group.add_argument(
        '--kmer-sizes', '--kmer_sizes', '-k',
        help='Manually specify the kmer-sizes used by SPAdes during assembly. Space separated odd integer values '
             'and less than 128 or "auto"',
        dest='kmer_sizes',
        default=['auto'],
        nargs="+",
    )

    assemble_group.add_argument(
        '--min-cov-long', '--min_cov_long',
        help='Automatically include Flye contigs with long read coverage greater than or equal to this. \n'
             'High long read coverage during assembly indicates that the overlap layout consensus algorithm \n'
             'is more likely to be correct.',
        dest='min_cov_long',
        default=20
    )

    assemble_group.add_argument(
        '--min-cov-short', '--min_cov_short',
        help='Automatically include Flye contigs with short read coverage less than or equal to this. \n'
             'Low coverage via short reads indicates that metaSPAdes will not be able to better assemble this contig.',
        dest='min_cov_short',
        default=3
    )

    assemble_group.add_argument(
        '--exclude-contig-cov', '--exclude_contig_cov',
        help='Automatically exclude Flye contigs with long read coverage less than or equal to this \n'
             'and less than or equal to `--exclude-contig-size`',
        dest='exclude_contig_cov',
        default=100
    )

    assemble_group.add_argument(
        '--exclude-contig-size', '--exclude_contig_size',
        help='Automatically exclude Flye contigs with length less than or equal to this \n'
             'and long read coverage less than or equal to `--exclude-contig-cov`',
        dest='exclude_contig_size',
        default=25000
    )

    assemble_group.add_argument(
        '--include-contig-size', '--include_contig_size',
        help='Automatically include Flye contigs with length less than or equal to this',
        dest='include_contig_size',
        default=100000
    )

    #~#~#~#~#~#~#~#~#~#~#~#~#~   sub-parsers   ~#~#~#~#~#~#~#~#~#~#~#~#~#
    ##########################  ~ ASSEMBLE ~  ###########################

    assemble_options = subparsers.add_parser('assemble',
                                              description='Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.',
                                              formatter_class=CustomHelpFormatter,
                                              parents=[qc_group, assemble_group, short_read_group, long_read_group, binning_group, base_group],
                                              epilog=
        '''
                                        ......:::::: ASSEMBLE ::::::......

        aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont

        ''')



    assemble_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['complete_assembly'],
    )

    ##########################  ~ RECOVER ~   ###########################

    recover_options = subparsers.add_parser('recover',
                                            description='The aviary binning pipeline',
                                            formatter_class=CustomHelpFormatter,
                                            parents=[qc_group, assemble_group, short_read_group, long_read_group, binning_group, annotation_group, base_group],
                                            epilog=
    '''
                                           ......:::::: RECOVER ::::::......
    
    aviary recover --assembly scaffolds.fasta -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont

    ''')

    recover_options.add_argument(
        '-a', '--assembly',
        help='Optional FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        nargs=1,
        required=False,
    )

    recover_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['recover_mags'],
    )

    recover_options.add_argument(
        '--perform-strain-analysis', '--perform_strain_analysis',
        help='Specify whether to use Lorikeet on recovered MAGs get strain diversity metrics',
        type=str2bool,
        nargs='?',
        const=True,
        dest='strain_analysis',
        default=False
    )

    ##########################  ~ ANNOTATE ~   ###########################

    annotate_options = subparsers.add_parser('annotate',
                                              description='Annotate a given set of MAGs using EggNOG, GTDB-tk, and Checkm2',
                                              formatter_class=CustomHelpFormatter,
                                              parents=[mag_group, annotation_group, base_group, qc_group],
                                              epilog=
                                            '''
                                                  ......:::::: ANNOTATE ::::::......
                                        
                                            aviary annotate --genome-fasta-directory input_bins/
                                        
                                            ''')

    annotate_options.add_argument(
        '-a', '--assembly',
        help='FASTA file containing scaffolded contigs of one or more metagenome assemblies wishing to be passed to QUAST',
        dest="assembly",
        nargs="*",
        required=False,
    )

    annotate_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['annotate'],
    )

    ##########################  ~ diversity ~   ###########################

    diversity_options = subparsers.add_parser('diversity',
                                             description='Perform strain diversity analysis',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[mag_group, qc_group, assemble_group, short_read_group, long_read_group,
                                                      binning_group, annotation_group, base_group],
                                             epilog=
                                             '''
                                                                    ......:::::: DIVERSITY ::::::......

                                             aviary diversity -c R1.fastq.gz R2.fastq.gz --genome-fasta-directory input_bins/

                                             ''')

    diversity_options.add_argument(
        '-a', '--assembly',
        help='FASTA file containing scaffolded contigs of one or more metagenome assemblies wishing to be passed to QUAST',
        dest="assembly",
        nargs="*",
        required=False,
    )

    diversity_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['lorikeet'],
    )

    diversity_options.add_argument(
        '--perform-strain-analysis', '--perform_strain_analysis',
        help=argparse.SUPPRESS,
        type=str2bool,
        nargs='?',
        const=True,
        dest='strain_analysis',
        default=True
    )

    ##########################  ~ CLUSTER ~   ###########################

    cluster_options = subparsers.add_parser('cluster',
                                             description='Clusters previous aviary runs together and performs'
                                                         'dereplication using Galah',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[base_group, cluster_group],
                                             epilog=
                                             '''
                                                                   ......:::::: CLUSTER ::::::......

                                             aviary cluster --input-runs aviary_output_folder_1/ aviary_output_folder_2/

                                             ''')

    cluster_options.add_argument(
        '-i', '--input-runs', '--input_runs',
        help='The paths to the previous finished runs of Aviary. Must contain the bins/checkm.out and bins/final_bins'
             'outputs',
        dest='previous_runs',
        nargs='*',
        required=True,
    )

    cluster_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['complete_cluster'],
    )

    ##########################  ~ VIRAL ~   ###########################

    viral_options = subparsers.add_parser('viral',
                                          description='The incomplete binning pipeline',
                                          formatter_class=CustomHelpFormatter,
                                          parents=[mag_group, short_read_group, long_read_group, annotation_group, base_group],
                                          epilog=
                                          '''
                                                  ......:::::: VIRAL ::::::...... 
 
                                          aviary viral --genome-fasta-files *.fasta
 
                                          ''')

    viral_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['create_webpage_genotype'],
    )

    ##########################   ~ COMPLETE ~  ###########################

    complete_options = subparsers.add_parser('complete',
                                            description='Performs all steps in the Aviary pipeline. '
                                                        'Assembly > Binning > Refinement > Annotation > Diversity',
                                            formatter_class=CustomHelpFormatter,
                                            parents=[qc_group, assemble_group, short_read_group, long_read_group, binning_group, annotation_group, base_group],
                                            epilog=
                                            '''
                                                               ......:::::: COMPLETE ::::::......

                                            aviary complete -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz 

                                            ''')

    complete_options.add_argument(
        '-a', '--assembly',
        help='Optional FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        nargs=1,
        required=False,
    )

    complete_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        nargs="+",
        default=['get_bam_indices', 'recover_mags', 'annotate', 'lorikeet'],
    )

    ##########################  ~ ISOLATE ~  ###########################

    isolate_options = subparsers.add_parser('isolate',
                                             description='Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[qc_group, short_read_group, long_read_group, isolate_group, binning_group, annotation_group, base_group],
                                             epilog=
                                             '''
                                                                             ......:::::: ISOLATE ::::::......
                                 
                                             aviary isolate -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont
                                 
                                             ''')

    isolate_options.add_argument(
        '-w', '--workflow',
        help='Main workflows to run',
        dest='workflow',
        nargs="+",
        default=['circlator'],
    )

    ##########################   ~ BATCH ~  ###########################

    batch_options = subparsers.add_parser('batch',
                                             description='Performs all steps in the Aviary pipeline on a batch file. \n'
                                                         'Each line in the batch file is processed separately and then \n'
                                                         'clustered using aviary. \n'
                                                         '(Assembly > Binning > Refinement > Annotation > Diversity) * n_samples --> Cluster',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[qc_group, assemble_group, binning_group, annotation_group, cluster_group, base_group],
                                             epilog=
                                             '''
                                                      ......:::::: BATCH ::::::......

                                             aviary batch -f batch_file.tsv -t 32 -o batch_test
                                             
                                             An example batch file can be found at: 

                                             ''')

    batch_options.add_argument(
        '-f', '--batch_file', '--batch-file',
        help='The tab or comma separated batch file containing the input samples to assemble and/or recover MAGs from. \n'
             'An example batch file can be found at XXX. The heading line is required. \n'
             'The number of reads provided to each sample is flexible as is the type of assembly being performed (if any). \n'
             'Multiple reads can be supplied by providing a comma-separated list (surrounded by double quotes \"\" if using a \n'
             'comma separated batch file) within the specific read column.',
        dest="batch_file",
        # nargs=1,
        required=True,
    )

    batch_options.add_argument(
        '--write-script', '--write_script',
        help='Write the aviary batch Snakemake commands to a bash script and exit. \n'
             'Useful when submitting jobs to HPC cluster with custom queueing.',
        dest='write_script',
        required=False
    )

    batch_options.add_argument(
        '--cluster',
        help='Cluster final output of all samples using aviary cluster if possible.',
        dest='cluster',
        type=str2bool,
        nargs='?',
        const=True,
        default=True
    )

    batch_options.add_argument(
        '--cluster-ani-values', '--cluster_ani_values', '--ani-values', '--ani_values',
        help='The range of ANI values to perform clustering and dereplication at during aviary cluster.',
        dest='ani_values',
        nargs='*',
        default=[0.99, 0.97, 0.95]
    )

    batch_options.add_argument(
        '--min-percent-read-identity-long', '--min_percent_read_identity_long',
        help='Minimum percent read identity used by CoverM for long-reads'
             'when calculating genome abundances.',
        dest='long_percent_identity',
        default='85'
    )

    batch_options.add_argument(
        '--min-percent-read-identity-short', '--min_percent_read_identity_short',
        help='Minimum percent read identity used by CoverM for short-reads \n'
             'when calculating genome abundances.',
        dest='short_percent_identity',
        default='95'
    )

    batch_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run for each sample',
        dest='workflow',
        nargs="+",
        default=['get_bam_indices', 'recover_mags', 'annotate', 'lorikeet'],
    )

    ##########################   ~ configure ~  ###########################

    configure_options = subparsers.add_parser('configure',
                                            description='Sets the conda environment variables for future runs and downloads databases. ',
                                            formatter_class=CustomHelpFormatter,
                                            parents=[base_group],
                                            epilog=
                                            '''
                                                               ......:::::: CONFIGURE ::::::......

                                            aviary configure --conda-prefix ~/.conda --gtdb-path ~/gtdbtk/release207/ --temp-dir /path/to/new/temp

                                            ''')

    configure_options.add_argument(
        '--gtdb-path', '--gtdb_path',
        help='Path to the local gtdb database files',
        dest='gtdb_path',
        required=False,
    )

    configure_options.add_argument(
        '--busco-db-path', '--busco_db_path',
        help='Path to the local BUSCO database files',
        dest='busco_db_path',
        required=False,
    )

    configure_options.add_argument(
        '--checkm2-db-path', '--checkm2_db_path',
        help='Path to Checkm2 Database',
        dest='checkm2_db_path',
        required=False,
    )

    configure_options.add_argument(
        '--eggnog-db-path', '--eggnog_db_path',
        help='Path to the local eggnog database files',
        dest='eggnog_db_path',
        required=False,
    )

    configure_options.add_argument(
        '-w', '--workflow',
        help=argparse.SUPPRESS,
        dest='workflow',
        nargs="+",
        default=['download_databases'],
    )

    ###########################################################################
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    args = main_parser.parse_args()
    time = datetime.now().strftime('%H:%M:%S %d-%m-%Y')

    if args.log:
        if os.path.isfile(args.log):
            raise Exception("File %s exists" % args.log)
        logging.basicConfig(filename=args.log,
                            level=debug[args.verbosity],
                            format='%(asctime)s %(levelname)s: %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=debug[args.verbosity],
                            format='%(asctime)s %(levelname)s: %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info("Time - %s" % (time))
    logging.info("Command - %s" % ' '.join(sys.argv))
    logging.info("Version - %s" % __version__)

    if args.subparser_name == 'configure':
        # Set the environment variables if manually configuring
        if args.conda_prefix is not None:
            Config.set_db_path(args.conda_prefix, db_name='CONDA_ENV_PATH')

        if args.tmpdir is not None:
            Config.set_db_path(args.tmpdir, db_name='TMPDIR')

        if args.gtdb_path is not None:
            Config.set_db_path(args.gtdb_path, db_name='GTDBTK_DATA_PATH')

        if args.busco_db_path is not None:
            Config.set_db_path(args.busco_db_path, db_name='BUSCO_DB')

        if args.checkm2_db_path is not None:
            Config.set_db_path(args.checkm2_db_path, db_name='CHECKM2DB')

        if args.eggnog_db_path is not None:
            Config.set_db_path(args.eggnog_db_path, db_name='EGGNOG_DATA_DIR')

        logging.info("The current aviary environment variables are:")
        logging.info(f"CONDA_ENV_PATH: {Config.get_software_db_path('CONDA_ENV_PATH', '--conda-prefix')}")
        logging.info(f"TMPDIR: {Config.get_software_db_path('TMPDIR', '--tmpdir')}")
        logging.info(f"GTDBTK_DATA_PATH: {Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')}")
        logging.info(f"EGGNOG_DATA_DIR: {Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')}")
        logging.info(f"CHECKM2DB: {Config.get_software_db_path('CHECKM2DB', '--checkm2-db-path')}")
        if not args.download:
            logging.info("All paths set. Exitting without downloading databases. If you wish to download databases use --download")
            sys.exit(0)

    # else:
    args = manage_env_vars(args)
    prefix = args.output
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    if args.subparser_name != 'batch':
        processor = Processor(args)
        processor.make_config()

        if args.build:
            try:
                args.cmds = args.cmds + '--conda-create-envs-only '
            except TypeError:
                args.cmds = '--conda-create-envs-only '

        try:
            if args.subparser_name == 'assemble':
                if args.use_unicycler:
                    args.workflow.insert(0, "combine_assemblies")
        except AttributeError:
            pass

        processor.run_workflow(cores=int(args.n_cores),
                               dryrun=args.dryrun,
                               clean=args.clean,
                               conda_frontend=args.conda_frontend,
                               snakemake_args=args.cmds)
    else:
        process_batch(args, prefix)

def manage_env_vars(args):
    if args.conda_prefix is None:
        args.conda_prefix = Config.get_software_db_path('CONDA_ENV_PATH', '--conda-prefix')

    if args.tmpdir is None:
        args.tmpdir = Config.get_software_db_path('TMPDIR', '--tmpdir')

    try:
        if args.gtdb_path is None:
            args.gtdb_path = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
        if args.eggnog_db_path is None:
            args.eggnog_db_path = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
        if args.checkm2_db_path is None:
            args.checkm2_db_path = Config.get_software_db_path('CHECKM2DB', '--checkm2-db-path')
    except AttributeError:
        pass

    return args

###############################################################################
################################ - Classes - ##################################

class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

if __name__ == '__main__':
    sys.exit(main())