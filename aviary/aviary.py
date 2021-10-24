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
from aviary.__init__ import __version__
import aviary.config.config as Config
from aviary.modules.processor import Processor
__author__ = "Rhys Newell"
__copyright__ = "Copyright 2020"
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
        cluster   - Clusters samples based on OTU content using SingleM **TBC**
        assemble  - Perform hybrid assembly using short and long reads, 
                    or assembly using only short reads
        recover   - Recover MAGs from provided assembly using a variety 
                    of binning algorithms 
        annotate  - Annotate MAGs **TBC**
        genotype  - Perform strain level analysis of MAGs **TBC**
        complete  - Runs each stage of the pipeline: assemble, recover, 
                    annotate, genotype in that order.

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
        help='The number of threads given to pplacer, values above 48 will be scaled down',
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
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default=Config.get_conda_path(),
    )

    base_group.add_argument(
        '--gtdb-path', '--gtdb_path',
        help='Path to the local gtdb files',
        dest='gtdb_path',
        default=Config.get_gtdb_path(),
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
        help='Which conda frontend to use, mamba is faster but harder to debug. Switch this to conda '
             'If experiencing problems installing environments',
        dest='conda_frontend',
        default="mamba",
        choices=["conda", "mamba"],
    )

    base_group.add_argument(
        '--clean',
        help='Clean up all temporary files. This will remove most BAM files and any FASTQ files '
             'generated from read filtering. Setting this to False is the equivalent of the --notemp'
             'option in snakemake. Useful for when running only part of a workflow as it avoids'
             'deleting files that would likely be needed in later parts of the workflow.'
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
        '--snakemake-cmds',
        help='Additional commands to supplied to snakemake in the form of a single string'
             'e.g. "--print-compilation True". '
             'NOTE: Most commands in snakemake -h are valid but some commands may clash with commands '
             'aviary directly supplies to snakemake. Please make'
             "sure your additional commands don't clash.",
        dest='cmds',
        default='',
    )

    ####################################################################
    qc_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                               add_help=False)
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
        default=1000
    )

    qc_group.add_argument(
        '--min-mean-q', '--min_mean_q',
        help='Minimum mean quality threshold',
        dest="min_mean_q",
        default=80
    )

    qc_group.add_argument(
        '--keep-percent', '--keep_percent',
        help='Percentage of reads passing quality thresholds kept by filtlong',
        dest="keep_percent",
        default=99
    )


    ####################################################################

    short_read_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                         add_help=False)
    read_group_exclusive = short_read_group.add_mutually_exclusive_group()

    read_group_exclusive.add_argument(
        '-1', '--pe-1', '--paired-reads-1', '--paired_reads_1', '--pe1',
        help='A space separated list of forwards read files to use for the binning process'
             'NOTE: If performing assembly and multiple files are provided then only '
             'the first file will be used for assembly.',
        dest='pe1',
        nargs='*',
        default="none"
    )

    short_read_group.add_argument(
        '-2', '--pe-2', '--paired-reads-2', '--paired_reads_2', '--pe2',
        help='A space separated list of forwards read files to use for the binning process'
             'NOTE: If performing assembly and multiple files are provided then only '
             'the first file will be used for assembly.',
        dest='pe2',
        nargs='*',
        default="none"
    )

    read_group_exclusive.add_argument(
        '-i','--interleaved',
        help='A space separated list of interleaved read files for the binning process '
             'NOTE: If performing assembly and multiple files are provided then only '
             'the first file will be used for assembly.',
        dest='interleaved',
        nargs='*',
        default="none"
    )

    read_group_exclusive.add_argument(
        '-c', '--coupled',
        help='Forward and reverse read files in a coupled space separated list. '
             'NOTE: If performing assembly and multiple files are provided then only '
             'the first two files will be used for assembly.',
        dest='coupled',
        nargs='*',
        default="none"
    )

    ####################################################################

    long_read_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                              add_help=False)
    long_read_group.add_argument(
        '-l', '--longreads', '--long-reads', '--long_reads',
        help='A space separated list of interleaved read files for the binning process. NOTE: If performing assembly and '
             'multiple long read files are provided, then only the first file is used for assembly. ',
        dest='longreads',
        nargs='*',
        default="none"
    )

    long_read_group.add_argument(
        '-z', '--longread-type', '--longread_type', '--long_read_type', '--long-read-type',
        help='Whether the sequencing platform and technology for the longreads. '
             '"rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS '
             'reads and "ont" for Oxford Nanopore',
        dest='longread_type',
        nargs=1,
        default="ont",
        choices=["ont", "rs", "sq", "ccs"],
    )

    ####################################################################

    annotation_group = argparse.ArgumentParser(add_help=False)

    annotation_group.add_argument(
        '--enrichm-db-path', '--enrichm_db_path',
        help='Path to the local EnrichM Database files',
        dest='enrichm_db_path',
        default=Config.get_enrichm_db_path(),
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

    ####################################################################
    mag_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                        add_help=False)
    mag_group_exclusive = mag_group.add_mutually_exclusive_group()

    mag_group_exclusive.add_argument(
        '-f', '--genome-fasta-files', '--genome_fasta_files',
        help='MAGs to be annotated',
        dest='mags',
        nargs='*',
        required=False,
    )

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
        default='fa'
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
        '-g', '--genome-size', '--genome_size',
        help='Approximate size of the isolate genome to be assembled',
        dest='genome_size',
        nargs=1,
        required=False,
        default=5000000
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
    #     '-g', '--genome-size', '--genome_size',
    #     help='Approximate size of the isolate genome to be assembled',
    #     dest='genome_size',
    #     nargs=1,
    #     required=False,
    #     default=5000000
    # )

    #~#~#~#~#~#~#~#~#~#~#~#~#~   sub-parsers   ~#~#~#~#~#~#~#~#~#~#~#~#~#
    ##########################   ~ CLUSTER ~  ###########################

    cluster_options = subparsers.add_parser('cluster',
                                             description='Cluster samples together based on OTU content. '
                                                         'Samples that cluster together should be used for assembly and binning.',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[short_read_group, long_read_group, base_group],
                                             epilog=
                                             '''
                                                                ......:::::: CLUSTER ::::::......
                                 
                                             aviary cluster -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz 
                                 
                                             ''')

    cluster_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='cluster_samples',
    )

    ##########################  ~ ASSEMBLE ~  ###########################

    assemble_options = subparsers.add_parser('assemble',
                                              description='Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.',
                                              formatter_class=CustomHelpFormatter,
                                              parents=[qc_group, short_read_group, long_read_group, binning_group, base_group],
                                              epilog=
        '''
                                        ......:::::: ASSEMBLE ::::::......

        aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont

        ''')


    assemble_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='complete_assembly',
    )

    ##########################  ~ RECOVER ~   ###########################

    recover_options = subparsers.add_parser('recover',
                                            description='The complete binning pipeline',
                                            formatter_class=CustomHelpFormatter,
                                            parents=[qc_group, short_read_group, long_read_group, binning_group, base_group],
                                            epilog=
    '''
                                           ......:::::: RECOVER ::::::......
    
    aviary recover --assembly scaffolds.fasta -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont

    ''')

    recover_options.add_argument(
        '-a', '--assembly',
        help='FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        nargs=1,
        required=False,
    )

    recover_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='recover_mags',
    )

    ##########################  ~ ANNOTATE ~   ###########################

    annotate_options = subparsers.add_parser('annotate',
                                              description='The complete binning pipeline',
                                              formatter_class=CustomHelpFormatter,
                                              parents=[mag_group, annotation_group, base_group],
                                              epilog=
                                            '''
                                                  ......:::::: ANNOTATE ::::::......
                                        
                                            aviary annotate --genome-fasta-files *.fasta
                                        
                                            ''')

    annotate_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='complete_annotation',
    )

    ##########################  ~ GENOTYPE ~   ###########################

    genotype_options = subparsers.add_parser('genotype',
                                             description='The complete binning pipeline',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[mag_group, short_read_group, long_read_group, base_group],
                                             epilog=
                                             '''
                                                     ......:::::: GENOTYPE ::::::......

                                             aviary genotype --genome-fasta-files *.fasta

                                             ''')

    genotype_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='create_webpage_genotype',
    )

    ##########################  ~ VIRAL ~   ###########################

    viral_options = subparsers.add_parser('viral',
                                              description='The complete binning pipeline',
                                              formatter_class=CustomHelpFormatter,
                                              parents=[mag_group, short_read_group, long_read_group, base_group],
                                              epilog=
                                              '''
                                                      ......:::::: VIRAL ::::::...... 
     
                                              aviary viral --genome-fasta-files *.fasta
     
                                              ''')

    viral_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='create_webpage_genotype',
    )

    ##########################   ~ COMPLETE ~  ###########################

    complete_options = subparsers.add_parser('complete',
                                            description='Cluster samples together based on OTU content. '
                                                        'Samples that cluster together should be used for assembly and binning.',
                                            formatter_class=CustomHelpFormatter,
                                            parents=[short_read_group, long_read_group, binning_group, annotation_group, base_group],
                                            epilog=
                                            '''
                                                               ......:::::: COMPLETE ::::::......

                                            aviary complete -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz 

                                            ''')

    complete_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='complete_workflow',
    )

    ##########################  ~ ISOLATE ~  ###########################

    isolate_options = subparsers.add_parser('isolate',
                                             description='Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.',
                                             formatter_class=CustomHelpFormatter,
                                             parents=[qc_group, short_read_group, long_read_group, isolate_group, binning_group, base_group],
                                             epilog=
                                             '''
                                                                             ......:::::: ISOLATE ::::::......
                                 
                                             aviary isolate -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont
                                 
                                             ''')

    isolate_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='create_webpage_assemble',
    )

    ##########################   ~ configure ~  ###########################

    configure_options = subparsers.add_parser('configure',
                                            description='Sets the conda environment variables for future runs. ',
                                            formatter_class=CustomHelpFormatter,
                                            epilog=
                                            '''
                                                               ......:::::: CONFIGURE ::::::......

                                            aviary configure --conda-prefix ~/.conda --gtdb-path ~/gtdbtk/release202/ 

                                            ''')

    configure_options.add_argument(
        '--conda-prefix', '--conda_prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        required=False,
    )

    configure_options.add_argument(
        '--gtdb-path', '--gtdb_path',
        help='Path to the local gtdb files',
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
        '--enrichm-db-path', '--enrichm_db_path',
        help='Path to the local EnrichM database files',
        dest='enrichm_db_path',
        required=False,
    )

    ###########################################################################
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    else:
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
                Config.set_conda_path(args.conda_prefix)

            if args.gtdb_path is not None:
                Config.set_gtdb_path(args.gtdb_path)

            if args.busco_db_path is not None:
                Config.set_busco_db_path(args.busco_db_path)

            if args.enrichm_db_path is not None:
                Config.set_enrichm_db_path(args.enrichm_db_path)

        else:
            prefix = args.output
            if not os.path.exists(prefix):
                os.makedirs(prefix)

            processor = Processor(args,
                               args.gtdb_path,
                               args.conda_prefix)

            processor.make_config()

            if args.build:
                try:
                    args.cmds = args.cmds + '--conda-create-envs-only '
                except TypeError:
                    args.cmds = '--conda-create-envs-only '

            processor.run_workflow(workflow=args.workflow,
                                   cores=int(args.n_cores),
                                   dryrun=args.dryrun,
                                   clean=args.clean,
                                   conda_frontend=args.conda_frontend,
                                   snakemake_args=args.cmds)


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
