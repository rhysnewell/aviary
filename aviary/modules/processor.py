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
import copy
from pathlib import Path
from glob import glob

# Local imports
from snakemake import utils
from snakemake.io import load_configfile
from ruamel.yaml import YAML  # used for yaml reading with comments
from aviary import LONG_READ_TYPES

BATCH_HEADER=['sample', 'short_reads_1', 'short_reads_2', 'long_reads', 'long_read_type', 'assembly', 'coassemble']

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
        self.tmpdir = os.path.abspath(args.tmpdir) if args.tmpdir else None
        self.resources = args.resources
        self.output = os.path.abspath(args.output)
        self.threads = args.max_threads
        self.max_memory = args.max_memory
        self.workflows = args.workflow
        self.request_gpu = args.request_gpu

        try:
            self.pplacer_threads = min(int(args.pplacer_threads), int(self.threads), 48)
        except AttributeError:
            self.pplacer_threads = min(int(self.threads), 48)

        try:
            self.strain_analysis = args.strain_analysis
        except AttributeError:
            self.strain_analysis = False

        # binning group items
        try:
            self.min_contig_size = args.min_contig_size
            self.min_bin_size = args.min_bin_size
            self.semibin_model = args.semibin_model
            self.refinery_max_iterations = args.refinery_max_iterations
            self.refinery_max_retries = args.refinery_max_retries
            self.skip_abundances = args.skip_abundances
            self.skip_taxonomy = args.skip_taxonomy
            self.skip_singlem = args.skip_singlem
            if args.binning_only:
                self.skip_abundances = True
                self.skip_taxonomy = True
                self.skip_singlem = True
            self.binning_only = args.binning_only

            self.skip_binners = ["maxbin2", "concoct", "comebin", "taxvamb"]
            if args.extra_binners:
                for binner in args.extra_binners:
                    binner = binner.lower()   
                    if binner == "maxbin" or binner == "maxbin2":
                        self.skip_binners.remove("maxbin2")
                    elif binner == "concoct":
                        self.skip_binners.remove("concoct")
                    elif binner == "comebin":
                        self.skip_binners.remove("comebin")
                    elif binner == "taxvamb":
                        self.skip_binners.remove("taxvamb")
                    else:
                        logging.warning(f"Unknown extra binner {binner} specified. Skipping...")

            if args.skip_binners:
                for binner in args.skip_binners:
                    binner = binner.lower()   
                    if binner == "metabat":
                        self.skip_binners.extend(["metabat_sens", "metabat_ssens", "metabat_spec", "metabat_sspec", "metabat2"])
                    elif binner == "metabat1":
                        self.skip_binners.extend(["metabat_sens", "metabat_ssens", "metabat_spec", "metabat_sspec"])
                    else:
                        self.skip_binners.append(binner)

        except AttributeError:
            self.min_contig_size = 1500
            self.min_bin_size = 200000
            self.semibin_model = 'global'
            self.refinery_max_iterations = 5
            self.refinery_max_retries = 3
            self.skip_binners = ["none"]
            self.skip_abundances = False
            self.binning_only = False
            self.skip_taxonomy = False
            self.skip_singlem = False

        try:
            self.assembly = args.assembly
        except AttributeError:
            self.assembly = 'none'

        try:
            self.reference_filter = [os.path.abspath(ref_fil) for ref_fil in args.reference_filter if ref_fil != 'none']
            if args.gold_standard is not None:
                self.gold_standard = [os.path.abspath(p) for p in args.gold_standard]
            else:
                self.gold_standard = 'none'
            
            self.min_read_size = args.min_read_size
            self.min_mean_q = args.min_mean_q
            self.keep_percent = args.keep_percent
            self.skip_qc = args.skip_qc
            self.min_short_read_size = args.min_short_read_length
            self.max_short_read_size = args.max_short_read_length
            self.disable_adapter_trimming = args.disable_adapter_trimming
            self.unqualified_percent_limit = args.unqualified_percent_limit
            self.quality_cutoff = args.quality_cutoff
            self.extra_fastp_params = args.extra_fastp_params
        except AttributeError:
            self.reference_filter = 'none'
            self.gold_standard = 'none'
            self.min_read_size = 0
            self.min_mean_q = 0
            self.keep_percent = 100
            self.skip_qc = False
            self.min_short_read_size = 0
            self.max_short_read_size = 0
            self.disable_adapter_trimming = False
            self.unqualified_percent_limit = 0
            self.quality_cutoff = 0
            self.extra_fastp_params = 'none'


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
            self.medaka_model = args.medaka_model
        except AttributeError:
            self.longread_type = 'none'
            self.medaka_model = 'none'

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
            self.coassemble = False
            self.min_cov_long = 20
            self.min_cov_short = 3
            self.exclude_contig_cov = 100
            self.exclude_contig_size = 25000
            self.long_contig_size = 100000

        try:
            self.mag_directory = os.path.abspath(args.directory) if args.directory is not None else 'none'
        except AttributeError:
            self.mag_directory = 'none'

        self.download = args.download

        try:
            if args.gtdb_path is not None:
                self.gtdbtk = args.gtdb_path
            else:
                self.gtdbtk = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
            if args.eggnog_db_path is not None:
                self.eggnog = args.eggnog_db_path
            else:
                self.eggnog = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
            if args.singlem_metapackage_path is not None:
                self.singlem = args.singlem_metapackage_path
            else:
                self.singlem = Config.get_software_db_path('SINGLEM_METAPACKAGE_PATH', '--singlem-metapackage-path')
            if args.metabuli_db_path is not None:
                self.metabuli = args.metabuli_db_path
            else:
                self.metabuli = Config.get_software_db_path('METABULI_DB_PATH', '--metabuli-db-path')
        except AttributeError:
            self.gtdbtk = Config.get_software_db_path('GTDBTK_DATA_PATH', '--gtdb-path')
            self.eggnog = Config.get_software_db_path('EGGNOG_DATA_DIR', '--eggnog-db-path')
            self.singlem = Config.get_software_db_path('SINGLEM_METAPACKAGE_PATH', '--singlem-metapackage-path')
            self.metabuli = Config.get_software_db_path('METABULI_DB_PATH', '--metabuli-db-path')
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
        
        if self.assembly == 'none' or self.assembly is None:
            # Check if coassembly or not needs to be specified by the user.
            if self.coassemble is None:
                if (self.pe1 != 'none' and len(self.pe1) > 1) or \
                   (self.longreads != 'none' and len(self.longreads) > 1):
                    logging.error("Multiple readsets detected. Either specify '--coassemble' for coassembly of or '--coassemble no'.")
                    sys.exit(-1)
        if self.coassemble is None:
            self.coassemble = False  # ensure that something is specified so that the config file is well formed

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
        conf["min_read_size"] = self.min_read_size
        conf["min_mean_q"] = self.min_mean_q
        conf["keep_percent"] = self.keep_percent
        conf["min_short_read_size"] = self.min_short_read_size
        conf["max_short_read_size"] = self.max_short_read_size
        conf["disable_adapter_trimming"] = self.disable_adapter_trimming
        conf["unqualified_percent_limit"] = self.unqualified_percent_limit
        conf["quality_cutoff"] = self.quality_cutoff
        conf["extra_fastp_params"] = self.extra_fastp_params
        conf["skip_qc"] = self.skip_qc
        conf["gsa"] = self.gold_standard
        conf["gsa_mappings"] = self.gsa_mappings
        conf["skip_binners"] = self.skip_binners
        conf["skip_abundances"] = self.skip_abundances
        conf["skip_taxonomy"] = self.skip_taxonomy
        conf["skip_singlem"] = self.skip_singlem
        conf["binning_only"] = self.binning_only
        conf["semibin_model"] = self.semibin_model
        conf["refinery_max_iterations"] = self.refinery_max_iterations
        conf["refinery_max_retries"] = self.refinery_max_retries
        conf["max_threads"] = int(self.threads)
        conf["pplacer_threads"] = int(self.pplacer_threads)
        conf["max_memory"] = int(self.max_memory)
        conf["request_gpu"] = self.request_gpu
        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.longreads
        conf["long_read_type"] = self.longread_type
        conf["medaka_model"] = self.medaka_model
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
        conf["download"] = self.download
        conf["gtdbtk_folder"] = self.gtdbtk
        conf["eggnog_folder"] = self.eggnog
        conf["singlem_metapackage"] = self.singlem
        conf["metabuli_folder"] = self.metabuli
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

    def run_workflow(self, cores=16, local_cores=None, profile=None, cluster_retries=None,
                     dryrun=False, clean=True,
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
        if self.tmpdir is not None:
            os.environ["TMPDIR"] = self.tmpdir
        for workflow in self.workflows:
            cmd = (
                "snakemake --snakefile {snakefile} --directory {working_dir} "
                "{jobs} {local_cores} --rerun-incomplete --keep-going {args} {rerun_triggers} "
                "--configfile {config_file} --nolock "
                "{profile} {retries} {conda_frontend} {resources} --use-conda {conda_prefix} "
                "{dryrun} {notemp} "
                "{target_rule}"
            ).format(
                snakefile=get_snakefile(),
                working_dir=self.output,
                jobs="--cores {}".format(cores) if cores is not None else "--jobs 1",
                local_cores="--local-cores {}".format(local_cores) if local_cores is not None else "",
                config_file=self.config,
                profile="" if not profile else "--profile {}".format(profile),
                retries="" if (cluster_retries is None) else "--retries {}".format(cluster_retries),
                dryrun="--dryrun" if dryrun else "",
                notemp="--notemp" if not clean else "",
                rerun_triggers="" if (rerun_triggers is None) else "--rerun-triggers {}".format(" ".join(rerun_triggers)),
                args=snakemake_args,
                target_rule=workflow if workflow != "None" else "",
                conda_prefix="--conda-prefix " + self.conda_prefix,
                conda_frontend="--conda-frontend " + "conda",
                resources=f"--resources mem_mb={int(self.max_memory)*1024} {self.resources}" if not dryrun else ""
            )

            logging.debug(f"Command: {cmd}")

            if write_to_script is not None:
                write_to_script.append(cmd)
                continue

            try:
                logging.info("Executing: %s" % cmd)
                subprocess.run(cmd.split(), check=True)
                logging.info("Finished: %s" % workflow)
                # logging.info("stderr: %s" % cmd_output)
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

    header=None
    separator=' '
    with open(args.batch_file, mode='r') as check_batch:
        for line in check_batch.readlines():
            line = line.strip()
            for sep in ['\t', ',', ' ']:
                separated = line.split(sep)
                if separated == BATCH_HEADER:
                    header=0
                    separator=sep
                    logging.debug("Inferred header")
                    break
                elif len(separated) >= 7:
                    header=None
                    separator=sep
                    logging.debug("Inferred no header")
                    break
            if header is None:
                logging.debug("No header found")
            break

    batch = pd.read_csv(args.batch_file, sep=separator, engine='python', names=BATCH_HEADER, header=header)
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
    except AttributeError:
        script_file = None
    
    write_to_script = None
    if script_file is not None:
        write_to_script = []

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
        if l_type not in LONG_READ_TYPES:
            logging.error(f"Unknown long read type {l_type} specified.")
            logging.error(f"Valid long read types: {LONG_READ_TYPES}")
            sys.exit(1)
        assembly = check_batch_input(batch.iloc[i, 5], None, split=False)
        coassemble = check_batch_input(batch.iloc[i, 6], False, split=False)
        
        new_args = copy.deepcopy(args)
        # update the value of args
        new_args.output = f"{prefix}/{sample}"
        runs.append(new_args.output)
        new_args.pe1 = s1
        new_args.pe2 = s2

        new_args.longreads = l
        new_args.longread_type = l_type
        new_args.assembly = assembly
        new_args.coassemble = coassemble

        # ensure output folder exists
        if not os.path.exists(new_args.output):
            os.makedirs(new_args.output)

        # setup processor for this line
        processor = Processor(new_args)
        processor.make_config()

        processor.run_workflow(cores=int(new_args.n_cores),
                               dryrun=new_args.dryrun,
                               clean=new_args.clean,
                               snakemake_args=new_args.cmds,
                               rerun_triggers=new_args.rerun_triggers,
                               profile=new_args.snakemake_profile,
                               cluster_retries=new_args.cluster_retries,
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
                                   snakemake_args=args.cmds,
                                   rerun_triggers=args.rerun_triggers,
                                   profile=args.snakemake_profile,
                                   cluster_retries=args.cluster_retries,
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
