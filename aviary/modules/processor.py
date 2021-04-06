#!/usr/bin/env python
###############################################################################
# processor.py - Class used to generate config file and make calls to the
#                snakemake pipeline
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
import logging
import os
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

###############################################################################
################################ - Functions - ################################

def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


# def update_config(config):
#     """
#     Populates config file with default config values.
#     And made changes if necessary.
#     """
#
#     # get default values and update them with values specified in config file
#     default_config = make_default_config()
#     utils.update_config(default_config, config)
#
#     return default_config


###############################################################################
################################ - Classes - ##################################

class Processor:
    def __init__(self,
                 args,
                 gtdbtk=Config.get_gtdb_path(),
                 conda_prefix=Config.get_conda_path(),
                 ):

        self.gtdbtk = gtdbtk
        self.conda_prefix = conda_prefix
        self.assembly = args.assembly
        self.reference_filter = args.reference_filter
        self.longreads = args.longreads
        self.longread_type = args.longread_type
        self.threads = args.max_threads
        self.max_memory = args.max_memory
        self.pplacer_threads = min(int(args.pplacer_threads), 48)
        self.min_contig_size = args.min_contig_size
        self.min_bin_size = args.min_bin_size
        self.output = args.output
        self.pe2 = args.pe2

        if args.interleaved == "none":
            self.pe1 = args.pe1
        elif args.pe2 == "none" and args.interleaved != "none":
            self.pe1 = args.interleaved

    def make_config(self):
        """
        Reads template config file with comments from ./template_config.yaml
        updates it by the parameters provided.
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
        conf["reference_filter"] = self.reference_filter
        conf["max_threads"] = self.threads
        conf["pplacer_threads"] = self.pplacer_threads
        conf["max_memory"] = int(self.max_memory)
        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type
        conf["min_contig_size"] = int(self.min_contig_size)
        conf["min_bin_size"] = int(self.min_bin_size)
        conf["gtdbtk_folder"] = os.path.abspath(self.gtdbtk)

        with open(self.config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
            "Configuration file written to %s\n" % self.config
        )

    def _validate_config(self):
        load_configfile(self.config)

    def run_workflow(self, workflow="recover_mags", cores=16, profile=None, dryrun=False, conda_frontend="mamba",
                     snakemake_args=""):
        """
        Runs the aviary pipeline
        By default all steps are executed
        Needs a config-file which is generated by given inputs.
        Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
        """

        if not os.path.exists(self.config):
            logging.critical(f"config-file not found: {self.config}\n")
            sys.exit(1)

        self._validate_config()

        cmd = (
            "snakemake --snakefile {snakefile} --directory {working_dir} "
            "{jobs} --rerun-incomplete "
            "--configfile '{config_file}' --nolock "
            " {profile} {conda_frontend} --use-conda {conda_prefix} {dryrun} "
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
            conda_frontend="--conda-frontend " + conda_frontend
        )
        logging.info("Executing: %s" % cmd)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            logging.critical(e)
            exit(1)