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

__author__ = "Rhys Newell"
__copyright__ = "Copyright 2020"
__credits__ = ["Rhys Newell"]
__license__ = "GPL3"
__version__ = "0.0.1"
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
        '--assembly',
        help='FASTA file containing scaffolded contigs of the metagenome assembly',
        dest="assembly",
        required=True,
    )

    input_options.add_argument(
        '--paired_reads_1',
        help='A space separated list of forwards read files to use for the binning process',
        dest='pe1',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '--paired_reads_2',
        help='A space separated list of forwards read files to use for the binning process',
        dest='pe2',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '--interleaved',
        help='A space separated list of interleaved read files for the binning process',
        dest='interleaved',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '--longreads',
        help='A space separated list of interleaved read files for the binning process',
        dest='longreads',
        nargs='*',
        default="none"
    )

    input_options.add_argument(
        '--longread_type',
        help='Whether the longreads are oxford nanopore or pacbio',
        dest='longread_type',
        nargs=1,
        default="nanopore",
        choices=["nanopore", "pacbio"],
    )

    input_options.add_argument(
        '--conda_prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default='~/.conda/envs/',
    )

    input_options.add_argument(
        '--gtdb_path',
        help='Path to the local gtdb files',
        dest='gtdb_path',
        default='/work/microbiome/db/gtdbtk/release95/',
    )

    input_options.add_argument(
        '--max_threads',
        help='Maximum number of threads given to any particular process',
        dest='max_threads',
        default=8,
    )

    input_options.add_argument(
        '--pplacer_threads',
        help='The number of threads given to pplacer, values above 48 will be scaled down',
        dest='pplacer_threads',
        default=8,
    )

    input_options.add_argument(
        '--n_cores',
        help='Maximum number of cores available for use. Must be >= to max_threads',
        dest='n_cores',
        default=16,
    )

    input_options.add_argument(
        '--output',
        help='Output directory, outputs to current directory *DON"T CHANGE, relative paths currently broken for this*',
        dest='output',
        default='./',
    )

    input_options.add_argument(
        '--workflow',
        help='Main workflow to run',
        dest='workflow',
        default='recover_mags',
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

        prefix = args.output
        if not os.path.exists(prefix):
            os.makedirs(prefix)
        if args.interleaved == "none":
            processor = aviary(args.assembly, args.pe1, args.pe2, args.longreads, args.longread_type,
                               int(args.max_threads), int(args.pplacer_threads), args.gtdb_path,
                               args.output, args.conda_prefix)
        elif args.pe2 == "none" and args.interleaved != "none":
            processor = aviary(args.assembly, args.interleaved, args.pe2, args.longreads, args.longread_type,
                               int(args.max_threads),
                               int(args.pplacer_threads),
                               args.gtdb_path,
                               args.output, args.conda_prefix)
        elif args.longreads != "none":
            processor = aviary(args.assembly, args.pe1, args.pe2,
                               args.longreads, args.longread_type,
                               int(args.max_threads),
                               int(args.pplacer_threads),
                               args.gtdb_path,
                               args.output, args.conda_prefix)
        else:
            sys.exit("Missing any input read files...")

        processor.make_config()
        processor.run_workflow(workflow=args.workflow, cores=int(args.n_cores))

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
                 gtdbtk="/work/microbiome/db/gtdbtk/release95",
                 output=".",
                 conda_prefix="~/.conda/envs/",
                 ):
        self.assembly = assembly
        self.pe1 = pe1
        self.pe2 = pe2
        self.longreads = longreads
        self.longread_type = longread_type
        self.threads = max_threads
        self.pplacer_threads = min(int(pplacer_threads), 48)
        self.gtdbtk = gtdbtk
        self.output = output
        self.conda_prefix = conda_prefix

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

        conf["fasta"] = self.assembly
        conf["max_threads"] = self.threads
        conf["pplacer_threads"] = self.pplacer_threads

        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type

        conf["gtdbtk_folder"] = self.gtdbtk


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
        For more details, see: https://metagenome-atlas.readthedocs.io
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
