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
import shutil
from datetime import datetime
import subprocess

# Local imports
from snakemake import utils
from snakemake.io import load_configfile
from ruamel.yaml import YAML  # used for yaml reading with comments

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

###############################################################################                                                                                                                      [44/1010]
################################ - Functions - ################################

def phelp():
    print(
    """
aviary

SUBCOMMAND:
recover
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

    ########################## ~ sub-parser ~ ###########################
    input_options = subparsers.add_parser('recover',
                                          description='The complete binning pipeline',
                                          formatter_class=CustomHelpFormatter,
                                          epilog='''
                                ~ RECOVER ~
    How to use recover:
    
    aviary recover --assembly scaffolds.fasta --paired_end_reads_1 *.1.fq.gz --paired_end_reads_1 *.2.fq.gz --longreads *.nanopore.fastq.gz --longread_type nanopore

    ''')

    input_options.add_argument(
        '-a', '--assembly',
        help='FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        nargs='*',
        required=True,
    )

    input_options.add_argument(
        '-1', '--paired-reads-1',
        help='A space separated list of forwards read files to use for the binning process',
        dest='pe1',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '-2','--paired-reads-2',
        help='A space separated list of forwards read files to use for the binning process',
        dest='pe2',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '-i','--interleaved',
        help='A space separated list of interleaved read files for the binning process',
        dest='interleaved',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '-l', '--longreads',
        help='A space separated list of interleaved read files for the binning process',
        dest='longreads',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '--longread-type',
        help='Whether the longreads are oxford nanopore or pacbio',
        dest='longread_type',
        nargs=1,
        default="nanopore",
        choices=["nanopore", "pacbio"],
    )

    input_options.add_argument(
        '-c', '--min-contig-size',
        help='Minimum contig size in base pairs to be considered for binning',
        dest='min_contig_size',
        nargs=1,
        default=1500
    )

    input_options.add_argument(
        '-s','--min-bin-size',
        help='Minimum bin size in base pairs for a MAG',
        dest='min_bin_size',
        nargs=1,
        default=200000
    )

    input_options.add_argument(
        '--conda-prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default='~/.conda/envs/',
    )

    input_options.add_argument(
        '--gtdb-path',
        help='Path to the local gtdb files',
        dest='gtdb_path',
        default='/work/microbiome/db/gtdbtk/release95/',
    )

    input_options.add_argument(
        '-t', '--max-threads',
        help='Maximum number of threads given to any particular process',
        dest='max_threads',
        default=8,
    )

    input_options.add_argument(
        '-p', '--pplacer-threads',
        help='The number of threads given to pplacer, values above 48 will be scaled down',
        dest='pplacer_threads',
        default=8,
    )

    input_options.add_argument(
        '-n', '--n-cores',
        help='Maximum number of cores available for use. Must be >= to max_threads',
        dest='n_cores',
        default=16,
    )

    input_options.add_argument(
        '-m', '--max-memory',
        help='Maximum memory for available usage in Gigabytes',
        dest='max_memory',
        default=250,
    )

    input_options.add_argument(
        '-o', '--output',
        help='Output directory',
        dest='output',
        default='./',
    )

    input_options.add_argument(
        '-w', '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='recover_mags',
    )

    input_options.add_argument(
        '--dry-run',
        help='Perform snakemake dry run, tests workflow order and conda environments',
        type=str2bool,
        nargs='?',
        const=True,
        dest='dryrun',
        default=False,
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

        prefix = args.output
        if not os.path.exists(prefix):
            os.makedirs(prefix)
        if args.interleaved == "none":
            processor = aviary(args.assembly,
                               args.pe1,
                               args.pe2,
                               args.longreads,
                               args.longread_type,
                               int(args.max_threads),
                               int(args.pplacer_threads),
                               int(args.max_memory),
                               args.gtdb_path,
                               args.output,
                               args.conda_prefix,
                               args)
        elif args.pe2 == "none" and args.interleaved != "none":
            processor = aviary(args.assembly, args.interleaved, args.pe2, args.longreads, args.longread_type,
                               int(args.max_threads),
                               int(args.pplacer_threads),
                               int(args.max_memory),
                               args.gtdb_path,
                               args.output, args.conda_prefix, args)
        elif args.longreads != "none":
            processor = aviary(args.assembly, args.pe1, args.pe2,
                               args.longreads, args.longread_type,
                               int(args.max_threads),
                               int(args.pplacer_threads),
                               int(args.max_memory),
                               args.gtdb_path,
                               args.output, args.conda_prefix, args)
        else:
            sys.exit("Missing any input read files...")

        processor.make_config()
        processor.run_workflow(workflow=args.workflow, cores=int(args.n_cores), dryrun=args.dryrun)

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


def update_config(config):
    """
    Populates config file with default config values.
    And made changes if necessary.
    """

    # get default values and update them with values specified in config file
    default_config = make_default_config()
    utils.update_config(default_config, config)

    return default_config

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

class aviary:
    def __init__(self,
                 assembly="none",
                 pe1="none",
                 pe2="none",
                 longreads="none",
                 longread_type="nanopore",
                 max_threads=16,
                 pplacer_threads=16,
                 max_memory=250,
                 gtdbtk="/work/microbiome/db/gtdbtk/release95",
                 output=".",
                 conda_prefix="~/.conda/envs/",
                 args=None
                 ):
        self.assembly = assembly
        self.pe1 = pe1
        self.pe2 = pe2
        self.longreads = longreads
        self.longread_type = longread_type
        self.threads = max_threads
        self.max_memory = max_memory
        self.pplacer_threads = min(int(pplacer_threads), 48)
        self.gtdbtk = gtdbtk
        self.output = output
        self.conda_prefix = conda_prefix
        if args is not None:
            self.min_contig_size = args.min_contig_size
            self.min_bin_size = args.min_bin_size

    def make_config(self):
        """
        Reads template config file with comments from ./template_config.yaml
        updates it by the parameters provided.
        Args:
            config (str): output file path for yaml
            database_dir (str): location of downloaded databases
            threads (int): number of threads per node to utilize
            assembler (str): either spades or megahit
            data_type (str): this is either metagenome or metatranscriptome
        """

        self.config = os.path.join(self.output, 'template_config.yaml')

        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False

        template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "template_config.yaml")

        with open(template_conf_file) as template_config:
            conf = yaml.load(template_config)

        if self.assembly != "none":
            self.assembly = [os.path.abspath(p) for p in self.assembly]
        if self.pe1 != "none":
            self.pe1 = [os.path.abspath(p) for p in self.pe1]
        if self.pe2 != "none":
            self.pe2 = [os.path.abspath(p) for p in self.pe2]
        if self.longreads != "none":
            self.longreads = [os.path.abspath(p) for p in self.longreads]
            
        conf["fasta"] = self.assembly
        conf["max_threads"] = self.threads
        conf["pplacer_threads"] = self.pplacer_threads
        conf["max_memory"] = self.max_memory


        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type
        conf["min_contig_size"] = self.min_contig_size
        conf["min_bin_size"] = self.min_bin_size

        conf["gtdbtk_folder"] = os.path.abspath(self.gtdbtk)


        with open(self.config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
            "Configuration file written to %s\n"
            "You may want to edit it using any text editor." % self.config
        )

    def validate_config(self):
        load_configfile(self.config)


    def run_workflow(self, workflow="recover_mags", cores=16, profile=None, dryrun=False, snakemake_args = ""):
        """Runs the aviary pipeline
        By default all steps are executed
        Needs a config-file which is generated by given inputs.
        Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
        """

        if not os.path.exists(self.config):
            logging.critical(f"config-file not found: {self.config}\n")
            sys.exit(1)

        self.validate_config()

        conf = load_configfile(self.config)

        cmd = (
            "snakemake --snakefile {snakefile} --directory {working_dir} "
            "{jobs} --rerun-incomplete "
            "--configfile '{config_file}' --nolock "
            " {profile} --use-conda {conda_prefix} {dryrun} "
            " {target_rule} "
            " {args} "
        ).format(
            snakefile=get_snakefile(),
            working_dir=self.output,
            jobs="--jobs {}".format(cores) if cores is not None else "",
            config_file=self.config,
            profile="" if (profile is None) else "--profile {}".format(profile),
            dryrun="--dryrun" if dryrun else "",
            args=" ".join(snakemake_args),
            target_rule=workflow if workflow != "None" else "",
            conda_prefix="--conda-prefix " + self.conda_prefix,
        )
        logging.info("Executing: %s" % cmd)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            logging.critical(e)
            exit(1)

if __name__ == '__main__':

    sys.exit(main())
