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
from pathlib import Path
from glob import glob

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
                 # conda_prefix=Config.get_software_db_path('CONDA_ENV_PATH', '--conda-prefix'),
                 ):


        self.conda_prefix = args.conda_prefix
        self.tmpdir = os.path.abspath(args.tmpdir)
        self.resources = args.resources
        self.output = os.path.abspath(args.output)
        self.threads = args.max_threads
        self.max_memory = args.max_memory
        self.pplacer_threads = min(int(self.threads), 48)
        self.workflows = args.workflow

        try:
            self.strain_analysis = args.strain_analysis
        except AttributeError:
            self.strain_analysis = False

        # binning group items
        try:
            self.min_contig_size = args.min_contig_size
            self.min_bin_size = args.min_bin_size
            self.semibin_model = args.semibin_model

            if args.skip_binners is not None:
                # skipping binners is tricky and requires that bam indices are produced
                # so we explicitly place the bam indices rule in the DAG before any binning occurs
                self.skip_binners = [binner.lower() for binner in args.skip_binners]
                self.check_binners_to_skip()
                if not any(rule in self.workflows for rule in ["get_bam_indices", "complete_assembly", "complete_assembly_with_qc", "complete_cluster"]):
                    self.workflows.insert(0, 'get_bam_indices')
                elif any(rule in self.workflows for rule in ["complete_assembly", "complete_assembly_with_qc"]):
                    indices = [idx for idx, rule in enumerate(self.workflows) if rule in ["complete_assembly", "complete_assembly_with_qc"]]
                    if len(indices) > 1:
                        logging.critical(f"Rules 'complete_assembly' and 'complete_assembly_with_qc' both found in DAG or one found multiple times.")
                        logging.critical(f"Please revise your --workflow parameter.")
                        sys.exit(1)
                    self.workflows.insert(indices[0] + 1, 'get_bam_indices')

        except AttributeError:
            self.min_contig_size = 1500
            self.min_bin_size = 200000
            self.semibin_model = 'global'
            self.skip_binners = []

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
            self.gsa_mappings = args.gsa_mappings
        except AttributeError:
            self.gsa_mappings = 'none'

        try:
            self.longreads = args.longreads
            self.long_percent_identity = args.long_percent_identity
        except AttributeError:
            self.longreads = 'none'
            self.long_percent_identity = 'none'

        try:
            self.longread_type = args.longread_type
        except AttributeError:
            self.longread_type = 'none'

        try:
            self.short_percent_identity = args.short_percent_identity

            if args.coupled != "none":
                self.pe1 = args.coupled[::2]
                self.pe2 = args.coupled[1::2]
                if len(self.pe1) != len(self.pe2):
                    logging.error(f"Number of forward reads != Number of reverse reads. Current forward: {len(self.pe1)} reverse: {len(self.pe2)}")
                    sys.exit(-1)
            else:
                self.pe2 = args.pe2
                if args.interleaved == "none":
                    self.pe1 = args.pe1
                elif args.pe2 == "none" and args.interleaved != "none":
                    self.pe1 = args.interleaved
        except AttributeError:
            self.pe1 = 'none'
            self.pe2 = 'none'
            self.short_percent_identity = 'none'


        try:
            self.kmer_sizes = args.kmer_sizes
            self.use_megahit = args.use_megahit
            self.coassemble = args.coassemble
            self.min_cov_long = args.min_cov_long
            self.min_cov_short = args.min_cov_short
            self.exclude_contig_cov = args.exclude_contig_cov
            self.exclude_contig_size = args.exclude_contig_size
            self.long_contig_size = args.include_contig_size
        except AttributeError:
            self.kmer_sizes = ['auto']
            self.use_megahit = False
            self.coassemble = True
            self.min_cov_long = 20
            self.min_cov_short = 3
            self.exclude_contig_cov = 100
            self.exclude_contig_size = 25000
            self.long_contig_size = 100000

        try:
            self.mag_directory = os.path.abspath(args.directory) if args.directory is not None else 'none'
        except AttributeError:
            self.mag_directory = 'none'

        try:
            if args.gtdb_path is not None:
                self.gtdbtk = args.gtdb_path
            else:
                self.gtdbtk = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
            if args.eggnog_db_path is not None:
                self.eggnog = args.eggnog_db_path
            else:
                self.eggnog = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
        except AttributeError:
            self.gtdbtk = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
            self.eggnog = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
            # self.enrichm = Config.get_software_db_path('ENRICHM_DB', '--enrichm-db-path')

        try:
            self.mag_extension = args.ext
        except AttributeError:
            self.mag_extension = 'none'

        try:
            self.previous_runs = [os.path.abspath(run) for run in args.previous_runs]
        except AttributeError:
            self.previous_runs = 'none'

        # Aviary cluster arguments
        try:
            if args.min_completeness == 'none':
                self.min_completeness = ' '
            else:
                self.min_completeness = f'--min-completeness {args.min_completeness}'

            if args.max_contamination == 'none':
                self.max_contamination = ' '
            else:
                self.max_contamination = f'--max-contamination {args.max_contamination}'

            self.precluster_ani = fraction_to_percent(args.precluster_ani)
            self.ani = fraction_to_percent(args.ani)
            self.precluster_method = args.precluster_method
            self.use_checkm2_scores = args.use_checkm2_scores
            self.pggb_params = args.pggb_params
        except AttributeError:
            self.min_completeness = 'none'
            self.max_contamination = 'none'
            self.ani = 'none'
            self.precluster_ani = 'none'
            self.precluster_method = 'none'
            self.use_checkm2_scores = False
            self.pggb_params = 'none'

        try:
            if args.checkm2_db_path is not None:
                self.checkm2_db = args.checkm2_db_path
            else:
                self.checkm2_db = Config.get_software_db_path('CHECKM2DB', '--checkm2-db-path')
        except AttributeError:
            self.checkm2_db = Config.get_software_db_path('CHECKM2DB', '--checkm2-db-path')
            # self.checkm2_db = 'none'

        # Must be always be first workflow
        if args.download:
            self.workflows.insert(0, 'download_databases')


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
            self.assembly = list(dict.fromkeys([os.path.abspath(p) for p in self.assembly]))
        elif self.assembly is None:
            self.assembly = 'none'
            logging.warning("No assembly provided, assembly will be created using available reads...")
        if self.pe1 != "none":
            self.pe1 = list(dict.fromkeys([os.path.abspath(p) for p in self.pe1]))
        if self.pe2 != "none":
            self.pe2 = list(dict.fromkeys([os.path.abspath(p) for p in self.pe2]))
        if self.longreads != "none":
            self.longreads = list(dict.fromkeys([os.path.abspath(p) for p in self.longreads]))
        if self.gsa_mappings != "none":
            self.gsa_mappings = os.path.abspath(self.gsa_mappings)

        conf["fasta"] = self.assembly
        conf["reference_filter"] = self.reference_filter
        conf["gsa"] = self.gold_standard
        conf["gsa_mappings"] = self.gsa_mappings
        conf["semibin_model"] = self.semibin_model
        conf["max_threads"] = int(self.threads)
        conf["pplacer_threads"] = int(self.pplacer_threads)
        conf["max_memory"] = int(self.max_memory)
        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type
        conf["kmer_sizes"] = self.kmer_sizes
        conf["use_megahit"] = self.use_megahit
        conf["coassemble"] = self.coassemble
        conf["min_cov_long"] = self.min_cov_long
        conf["min_cov_short"] = self.min_cov_short
        conf["exclude_contig_cov"] = self.exclude_contig_cov
        conf["exclude_contig_size"] = self.exclude_contig_size
        conf["long_contig_size"] = self.long_contig_size
        conf["min_contig_size"] = int(self.min_contig_size)
        conf["min_bin_size"] = int(self.min_bin_size)
        conf["gtdbtk_folder"] = self.gtdbtk
        conf["eggnog_folder"] = self.eggnog
        conf["strain_analysis"] = self.strain_analysis
        conf["checkm2_db_folder"] = self.checkm2_db
        conf["use_checkm2_scores"] = self.use_checkm2_scores
        conf["mag_directory"] = self.mag_directory
        conf["mag_extension"] = self.mag_extension
        conf["previous_runs"] = self.previous_runs
        conf["min_completeness"] = self.min_completeness
        conf["max_contamination"] = self.max_contamination
        conf["ani"] = self.ani
        conf["precluster_ani"] = self.precluster_ani
        conf["precluster_method"] = self.precluster_method
        conf["pggb_params"] = self.pggb_params
        conf["tmpdir"] = self.tmpdir

        with open(self.config, "w") as f:
            yaml.dump(conf, f)
        logging.info(
            "Configuration file written to %s" % self.config
        )

    def _validate_config(self):
        load_configfile(self.config)

    def check_binners_to_skip(self):

        to_skip = []
        skipped = 0

        if "rosella" in self.skip_binners:
            rosella_path = self.output + "/data/rosella_bins/"
            rosella_ref_path = self.output + "/data/rosella_refined/"
            os.makedirs(rosella_path, exist_ok=True)
            os.makedirs(rosella_ref_path, exist_ok=True)
            to_skip.append(rosella_path + "done")
            to_skip.append(rosella_ref_path + "done")
            skipped += 1

        if "semibin" in self.skip_binners:
            sb_path = self.output + "/data/semibin_bins/"
            sb_ref_path = self.output + "/data/semibin_refined/"
            os.makedirs(sb_path, exist_ok=True)
            os.makedirs(sb_ref_path, exist_ok=True)
            to_skip.append(sb_path + "done")
            to_skip.append(sb_ref_path + "done")
            skipped += 1

        if "vamb" in self.skip_binners:
            vamb_path = self.output + "/data/vamb_bins/"
            os.makedirs(vamb_path, exist_ok=True)
            to_skip.append(vamb_path + "done")
            skipped += 1

        if "concoct" in self.skip_binners:
            concoct_path = self.output + "/data/concoct_bins/"
            os.makedirs(concoct_path, exist_ok=True)
            to_skip.append(concoct_path + "done")
            skipped += 1

        if any(m in self.skip_binners for m in ["maxbin", "maxbin2"]):
            maxbin2_path = self.output + "/data/maxbin2_bins/"
            os.makedirs(maxbin2_path, exist_ok=True)
            to_skip.append(maxbin2_path + "done")
            skipped += 1

        if any(m in self.skip_binners for m in ["metabat", "metabat1"]):
            to_skip.extend(self.skip_metabat1())
            skipped += 1

        if any(m in self.skip_binners for m in ["metabat", "metabat2"]):
            m_path = self.output + "/data/metabat_bins_2/"
            m_ref_path = self.output + "/data/metabat2_refined/"
            os.makedirs(m_path, exist_ok=True)
            os.makedirs(m_ref_path, exist_ok=True)
            to_skip.append(m_path + "done")
            to_skip.append(m_ref_path + "done")
            skipped += 1

        if skipped < 7:
            [Path(skip).touch(exist_ok=True) for skip in to_skip]
        else:
            logging.error("Check --skip-binners. All binners are being skipped. At least one binning algorithm must be used.")
            sys.exit(1)


    def skip_metabat1(self):
        return_paths = []
        for m in ["metabat_bins_sens/", "metabat_bins_ssens/", "metabat_bins_spec/", "metabat_bins_sspec/"]:
            m_path = self.output + "/data/" + m
            os.makedirs(m_path, exist_ok=True)
            return_paths.append(m_path + "done")

        return return_paths

    def run_workflow(self, cores=16, profile=None,
                     dryrun=False, clean=True, conda_frontend="mamba",
                     snakemake_args="", write_to_script=None, rerun_triggers=None):
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

        for workflow in self.workflows:
            cmd = (
                "snakemake --snakefile {snakefile} --directory {working_dir} "
                "{jobs} --rerun-incomplete {args} {rerun_triggers} "
                "--configfile '{config_file}' --nolock "
                "{profile} {conda_frontend} {resources} --use-conda {conda_prefix} "
                "{dryrun} {notemp} "
                "{target_rule}"
            ).format(
                snakefile=get_snakefile(),
                working_dir=self.output,
                jobs="--jobs {}".format(cores) if cores is not None else "",
                config_file=self.config,
                profile="" if (profile is None) else "--profile {}".format(profile),
                dryrun="--dryrun" if dryrun else "",
                notemp="--notemp" if not clean else "",
                rerun_triggers="" if (rerun_triggers is None) else "--rerun-triggers {}".format(" ".join(rerun_triggers)),
                args=snakemake_args,
                target_rule=workflow if workflow != "None" else "",
                conda_prefix="--conda-prefix " + self.conda_prefix,
                conda_frontend="--conda-frontend " + conda_frontend,
                resources=f"--default-resources \"tmpdir='{self.tmpdir}'\" --resources mem_mb={int(self.max_memory)*1024} {self.resources}" if not dryrun else ""
            )

            if write_to_script is not None:
                write_to_script.append(cmd)
                continue

            try:
                subprocess.check_call(cmd, shell=True)
                logging.info("Executing: %s" % cmd)
            except subprocess.CalledProcessError as e:
                # removes the traceback
                logging.critical(e)
                exit(1)


def process_batch(args, prefix):
    '''
    Function for handling and processing aviary commands in batches
    The user supplies a tab separated batch file with six defined columns. Aviary is run on each line
    of the batch file using the specified workflow.

    If the user wishes, the results are then clustered.
    '''
    import pandas as pd

    logging.info(f"Reading batch file: {args.batch_file}")

    header=0
    with open(args.batch_file, mode='r') as check_batch:
        for line in check_batch.readlines():
            if "sample\tshort_reads_1\tshort_reads_2\tlong_reads\tlong_read_type\tassembly\tcoassemble" in line \
                or "sample,short_reads_1,short_reads_2,long_reads,long_read_type,assembly,coassemble" in line \
                or "sample  short_reads_1   short_reads_2   long_reads      long_read_type  assembly        coassemble" in line \
                or "sample short_reads_1 short_reads_2 long_reads long_read_type assembly coassemble" in line:
               header=1
               logging.debug("Inferred header")
            else:
                logging.debug("No heading inferred.")
            break

    batch = pd.read_csv(args.batch_file, sep=None, engine='python', skiprows=header)
    if len(batch.columns) != 7:
        logging.critical(f"Batch file contains incorrect number of columns ({len(batch.columns)}). Should contain 7.")
        logging.critical(f"Current columns: {batch.columns}")
        sys.exit()

    if args.build:
        try:
            args.cmds = args.cmds + '--conda-create-envs-only '
        except TypeError:
            args.cmds = '--conda-create-envs-only '

    if args.use_unicycler:
        args.workflow.insert(0, "combine_assemblies")

    try:
        script_file = args.write_script
        write_to_script = []
    except AttributeError:
        script_file = None
        write_to_script = None

    runs = []
    args.interleaved = "none" # hacky solution to skip attribute error
    args.coupled = "none"
    for i in range(batch.shape[0]):
        # process the batch line
        sample = check_batch_input(batch.iloc[i, 0], f"sample_{i}", split=False)
        logging.info(f"Processing {sample}")
        s1 = check_batch_input(batch.iloc[i, 1], "none", split=True)
        s2 = check_batch_input(batch.iloc[i, 2], "none", split=True)
        l = check_batch_input(batch.iloc[i, 3], "none", split=True)
        l_type = check_batch_input(batch.iloc[i, 4], "ont", split=False)
        assembly = check_batch_input(batch.iloc[i, 5], None, split=False)
        coassemble = check_batch_input(batch.iloc[i, 6], False, split=False)

        # update the value of args
        args.output = f"{prefix}/{sample}"
        runs.append(args.output)
        args.pe1 = s1
        args.pe2 = s2

        args.longreads = l
        args.longread_type = l_type
        args.assembly = assembly
        args.coassemble = coassemble

        # ensure output folder exists
        if not os.path.exists(args.output):
            os.makedirs(args.output)

        # setup processor for this line
        processor = Processor(args)
        processor.make_config()

        processor.run_workflow(cores=int(args.n_cores),
                               dryrun=args.dryrun,
                               clean=args.clean,
                               conda_frontend=args.conda_frontend,
                               snakemake_args=args.cmds,
                               write_to_script=write_to_script)

    if args.cluster:
        logging.info(f"Beginning clustering of {len(runs)} previous Aviary runs with ANI values: {args.ani_values}...")

        for ani in args.ani_values:
            args.previous_runs = runs
            args.ani = ani
            args.workflow = ['complete_cluster']
            args.output = f"{prefix}/aviary_cluster_ani_{ani}"
            # ensure output folder exists
            if not os.path.exists(args.output):
                os.makedirs(args.output)
            processor = Processor(args)
            processor.make_config()

            processor.run_workflow(cores=int(args.n_cores),
                                   dryrun=args.dryrun,
                                   clean=args.clean,
                                   conda_frontend=args.conda_frontend,
                                   snakemake_args=args.cmds,
                                   write_to_script=write_to_script)

    if script_file is not None:
        with open(script_file, 'w') as sf:
            for line in write_to_script:
                sf.write(f"{line}\n")


def check_batch_input(val, default=None, split=False, split_val=','):
    """
    Takes a batch entry from within a line and ensures the output makes sense for Aviary
    Split is required for cells that separated by the split_val
    """
    if not isinstance(val, str):
        return default

    new_val = val.strip()

    if split:
        new_val = new_val.split(split_val)
        return_vals = []
        for paths in new_val:
            return_vals.extend(glob(paths))
        new_val = return_vals

    return new_val

def fraction_to_percent(val):
    val = float(val)
    if val <= 1:
        return val * 100
    return val