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
import glob

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
                 conda_prefix=Config.get_software_db_path('CONDA_ENV_PATH', '--conda-prefix'),
                 ):


        self.conda_prefix = conda_prefix

        try:
            self.assembly = args.assembly
        except AttributeError:
            self.assembly = 'none'

        try:

            self.reference_filter = os.path.abspath(args.reference_filter) if args.reference_filter != 'none' else 'none'
            if args.gold_standard is not None:
                self.gold_standard = [os.path.abspath(p) for p in args.gold_standard]
            else:
                self.gold_standard = 'none'
        except AttributeError:
            self.reference_filter = 'none'
            self.gold_standard = 'none'

        try:
            self.longreads = args.longreads
        except AttributeError:
            self.longreads = 'none'
        try:
            self.longread_type = args.longread_type
        except AttributeError:
            self.longread_type = 'none'
        self.threads = args.max_threads
        self.max_memory = args.max_memory
        self.pplacer_threads = min(int(self.threads), 48)

        try:
            self.min_contig_size = args.min_contig_size
        except AttributeError:
            self.min_contig_size = 1500

        try:
            self.min_bin_size = args.min_bin_size
        except AttributeError:
            self.min_bin_size = 200000
        self.output = args.output

        try:
            if args.coupled != "none":
                self.pe1 = args.coupled[::2]
                self.pe2 = args.coupled[1::2]
            else:
                self.pe2 = args.pe2
                if args.interleaved == "none":
                    self.pe1 = args.pe1
                elif args.pe2 == "none" and args.interleaved != "none":
                    self.pe1 = args.interleaved
        except AttributeError:
            self.pe1 = 'none'
            self.pe2 = 'none'

        try:
            self.mag_directory = os.path.abspath(args.directory) if args.directory is not None else 'none'
        except AttributeError:
            self.mag_directory = 'none'


        try:
            self.gtdbtk = args.gtdb_path
            self.eggnog = args.eggnog_db_path
            # self.enrichm = args.enrichm_db_path
        except AttributeError:
            self.gtdbtk = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
            self.eggnog = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
            # self.enrichm = Config.get_software_db_path('ENRICHM_DB', '--enrichm-db-path')
        # try:
        #     self.mags = args.mags
        # except AttributeError:
        #     self.mags = 'none'

        try:
            self.mag_extension = args.ext
        except AttributeError:
            self.mag_extension = 'none'

        try:
            self.previous_runs = [os.path.abspath(run) for run in args.previous_runs]
        except AttributeError:
            self.previous_runs = 'none'
        
        try:
            self.min_completeness = args.min_completeness
            self.max_contamination = args.max_contamination
            self.ani = args.ani
            self.precluster_ani = args.precluster_ani
            self.precluster_method = args.precluster_method
        except AttributeError:
            self.min_completeness = 'none'
            self.max_contamination = 'none'
            self.ani = 'none'
            self.precluster_ani = 'none'
            self.precluster_method = 'none'


        try:
            # self.checkm2_db = Config.get_software_db_path('CHECKM2DB', '--checkm2-db-path')
            self.checkm2_db = os.environ['CHECKM2DB']
        except KeyError:
            self.checkm2_db = 'none'

    def make_config(self):
        """
        Reads template config file with comments from ./template_config.yaml
        updates it by the parameters provided.
        """

        self.config = os.path.join(self.output, 'config.yaml')

        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False

        template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "template_config.yaml")

        with open(template_conf_file) as template_config:
            conf = yaml.load(template_config)

        if self.assembly != "none" and self.assembly is not None:
            self.assembly = [os.path.abspath(p) for p in self.assembly]
        elif self.assembly is None:
            self.assembly = 'none'
            logging.warning("No assembly provided, assembly will be created using available reads...")
        if self.pe1 != "none":
            self.pe1 = [os.path.abspath(p) for p in self.pe1]
        if self.pe2 != "none":
            self.pe2 = [os.path.abspath(p) for p in self.pe2]
        if self.longreads != "none":
            self.longreads = [os.path.abspath(p) for p in self.longreads]

        conf["fasta"] = self.assembly
        conf["reference_filter"] = self.reference_filter
        conf["gsa"] = self.gold_standard
        conf["max_threads"] = int(self.threads)
        conf["pplacer_threads"] = int(self.pplacer_threads)
        conf["max_memory"] = int(self.max_memory)
        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type
        conf["min_contig_size"] = int(self.min_contig_size)
        conf["min_bin_size"] = int(self.min_bin_size)
        conf["gtdbtk_folder"] = self.gtdbtk
        conf["eggnog_folder"] = self.eggnog
        # conf["enrichm_folder"] = self.enrichm

        conf["checkm2_db_path"] = self.checkm2_db
        conf["mag_directory"] = self.mag_directory
        conf["mag_extension"] = self.mag_extension
        # conf["mags"] = self.mags
        conf["previous_runs"] = self.previous_runs
        conf["min_completeness"] = self.min_completeness
        conf["max_contamination"] = self.max_contamination
        conf["ani"] = self.ani
        conf["precluster_ani"] = self.precluster_ani
        conf["precluster_method"] = self.precluster_method

        with open(self.config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
            "Configuration file written to %s\n" % self.config
        )

    def _validate_config(self):
        load_configfile(self.config)

    def run_workflow(self, workflow="recover_mags", cores=16, profile=None,
                     dryrun=False, clean=True, conda_frontend="mamba",
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

        cores = max(int(self.threads), cores)

        cmd = (
            "snakemake --snakefile {snakefile} --directory {working_dir} "
            "{jobs}--rerun-incomplete "
            "--configfile '{config_file}' --nolock"
            " {profile} {conda_frontend} --use-conda {conda_prefix}"
            " {dryrun}{notemp}{args}"
            " {target_rule}"
        ).format(
            snakefile=get_snakefile(),
            working_dir=self.output,
            jobs="--jobs {} ".format(cores) if cores is not None else "",
            config_file=self.config,
            profile="" if (profile is None) else "--profile {}".format(profile),
            dryrun="--dryrun " if dryrun else "",
            notemp="--notemp " if not clean else "",
            args=snakemake_args,
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