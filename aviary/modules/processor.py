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
import copy
from pathlib import Path
from glob import glob

# Local imports
from snakemake import utils
from snakemake.io import load_configfile
from ruamel.yaml import YAML  # used for yaml reading with comments
from aviary import LONG_READ_TYPES, COVERAGE_JOB_STRATEGIES, COVERAGE_JOB_CUTOFF
from aviary.modules.common import workflow_identifier
import re
import threading

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
    def __init__(self, args):
        self.tmpdir = os.path.abspath(args.tmpdir) if args.tmpdir else None
        self.resources = args.resources
        self.output = os.path.abspath(args.output)
        self.threads = args.max_threads
        self.max_memory = args.max_memory
        self.workflows = args.workflow
        self.request_gpu = args.request_gpu
        self.strict = args.strict

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

            if args.coverage_job_strategy == COVERAGE_JOB_STRATEGIES[0]:
                num_samples = 0
                if args.pe1 != "none":
                    num_samples += len(args.pe1)
                if args.interleaved != "none":
                    num_samples += len(args.interleaved)
                if args.longreads != "none":
                    num_samples += len(args.longreads)

                self.coverage_split = num_samples >= COVERAGE_JOB_CUTOFF
            elif args.coverage_job_strategy == COVERAGE_JOB_STRATEGIES[1]:
                self.coverage_split = True
            else:
                self.coverage_split = False

            self.coverage_samples_per_job = args.coverage_samples_per_job
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
            self.coverage_split = False
            self.coverage_samples_per_job = 5
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
            if args.host_filter != ['none']:
                self.host_filter = [os.path.abspath(ref_fil) for ref_fil in args.host_filter if ref_fil != 'none']
            else:
                self.host_filter = ['none']

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
            self.host_filter = ['none']
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

        # Ensure that all input read files we have read permission on
        if self.pe1 != 'none':
            for p in self.pe1:
                if not os.access(p, os.R_OK):
                    logging.error(f"Cannot read short read file {p}. Please check permissions.")
                    sys.exit(1)
        if self.pe2 != 'none':
            for p in self.pe2:
                if not os.access(p, os.R_OK):
                    logging.error(f"Cannot read short read file {p}. Please check permissions.")
                    sys.exit(1)
        if self.longreads != 'none':
            for p in self.longreads:
                if not os.access(p, os.R_OK):
                    logging.error(f"Cannot read long read file {p}. Please check permissions.")
                    sys.exit(1)

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
        conf["host_filter"] = self.host_filter
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
        conf["coverage_split"] = self.coverage_split
        conf["coverage_samples_per_split"] = self.coverage_samples_per_job
        conf["refinery_max_iterations"] = self.refinery_max_iterations
        conf["refinery_max_retries"] = self.refinery_max_retries
        conf["max_threads"] = int(self.threads)
        conf["pplacer_threads"] = int(self.pplacer_threads)
        conf["max_memory"] = int(self.max_memory)
        conf["request_gpu"] = self.request_gpu
        conf["strict"] = self.strict
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

    def parse_snakemake_errors(self, text: str):
        """Extract failed rule names from Snakemake output.
        Since Snakemake no longer reliably prints "log:" lines, we only
        capture the erroring rule names here; log discovery happens via the
        workflow's resources:log_path layout on disk.
        """
        rules = []
        for raw in text.splitlines():
            m = re.search(r"Error in rule (\S+):", raw)
            if m:
                rules.append(m.group(1))
        # Preserve order but dedupe
        seen = set()
        unique_rules = []
        for r in rules:
            if r not in seen:
                unique_rules.append(r)
                seen.add(r)
        return unique_rules

    def logs_dir_root_for_workflow(self, output_dir: str, workflow: str) -> str:
        # All module Snakefiles use logs_dir = "logs" relative to the working directory
        # (self.output). Therefore logs live in a single global folder: <output>/logs.
        return os.path.join(output_dir, "logs")

    def find_logs_for_rule(self, rule_name: str, workflow: str, output_dir: str, wid: str | None = None):
        """Discover log files for a rule by scanning module Snakefiles for its resources:log_path.
        Parameters
        - rule_name: Snakemake rule name
        - workflow: Target rule invoked in the global Snakefile (e.g., 'coassemble.smk')
        - output_dir: Snakemake working directory used by run_workflow
        - wid: workflow identifier subdirectory; defaults to aviary.modules.common.workflow_identifier
        """

        wid = wid or workflow_identifier
        logs_root = self.logs_dir_root_for_workflow(output_dir, workflow)
        if not os.path.isdir(logs_root):
            return []

        # Scan all module Snakefiles for the rule definition
        modules_dir = os.path.dirname(os.path.abspath(__file__))
        candidate_files = []
        # Include global Snakefile just in case, then all module *.smk
        candidate_files.append(os.path.join(modules_dir, "Snakefile"))
        for root, _, files in os.walk(modules_dir):
            for fn in files:
                if fn.endswith(".smk"):
                    candidate_files.append(os.path.join(root, fn))

        patterns = []
        rule_block_re = re.compile(rf"(?ms)^rule\s+{re.escape(rule_name)}\s*:\s*(.*?)\n(?=rule\s+\w+:|$)")
        log_lambda_re = re.compile(r"log_path\s*=\s*lambda[^:]*:\s*setup_log\((.+?),\s*attempt\)", re.S)
        fstring_re = re.compile(r"f[\"'](.+?)[\"']", re.S)

        for path in candidate_files:
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as sf:
                    text = sf.read()
            except Exception:
                continue

            m_block = rule_block_re.search(text)
            if not m_block:
                continue

            block = m_block.group(1)
            m_log = log_lambda_re.search(block)
            if not m_log:
                # Rule found but no log_path resource; keep scanning others, but we'll have fallback
                continue

            arg = m_log.group(1)
            m_f = fstring_re.search(arg)
            if not m_f:
                continue

            content = m_f.group(1)
            # Replace placeholders
            # {logs_dir} -> logs_root
            content = content.replace("{logs_dir}", logs_root)
            # Replace wildcards with glob star
            content = re.sub(r"\{wildcards\.[^}]+\}", "*", content)
            # Normalize duplicate slashes
            content = re.sub(r"/+", "/", content)
            # Build glob pattern to attempts; prefer any wid
            patterns.append(os.path.join(content, "*", "attempt*.log"))

            # Stop after first successful rule match; most specific
            break

        # Fallback if no specific pattern found for the rule
        if not patterns:
            patterns.append(os.path.join(logs_root, "**", "attempt*.log"))

        files = set()
        for pat in patterns:
            for f in glob(pat, recursive=True):
                if os.path.isfile(f):
                    files.add(os.path.abspath(f))

        # Prefer higher attempt numbers and newer wid directories. Keep deterministic order.
        def attempt_num(p):
            m = re.search(r"attempt(\d+)\.log$", os.path.basename(p))
            return int(m.group(1)) if m else -1
        def wid_dir_key(p):
            # parent dir name is the wid (e.g., 20250101_123456); lexical sort works
            return os.path.basename(os.path.dirname(p))

        return sorted(files, key=lambda p: (wid_dir_key(p), attempt_num(p)), reverse=True)

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
                "{jobs} {local_cores} --rerun-incomplete --keep-going {args} {rerun_triggers} {resources} "
                "--configfile {config_file} --nolock "
                "{profile} {retries} "
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
                resources=f"--resources mem_mb={int(self.max_memory)*1024} {self.resources}" if not dryrun else ""
            )

            logging.debug(f"Command: {cmd}")

            if write_to_script is not None:
                write_to_script.append(cmd)
                continue

            # Stream stdout and stderr separately in real time while buffering for parsing
            combined_output = []
            out_buffer = []
            err_buffer = []
            logging.info("Executing: %s" % cmd)
            proc = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
            )

            def _pump(stream, writer, buf):
                try:
                    for line in stream:
                        writer.write(line)
                        writer.flush()
                        buf.append(line)
                        combined_output.append(line)
                finally:
                    try:
                        stream.close()
                    except Exception:
                        pass

            threads = []
            if proc.stdout is not None:
                t_out = threading.Thread(target=_pump, args=(proc.stdout, sys.stdout, out_buffer))
                t_out.daemon = True
                t_out.start()
                threads.append(t_out)
            if proc.stderr is not None:
                t_err = threading.Thread(target=_pump, args=(proc.stderr, sys.stderr, err_buffer))
                t_err.daemon = True
                t_err.start()
                threads.append(t_err)

            for t in threads:
                t.join()
            proc.wait()

            if proc.returncode == 0:
                logging.info("Finished: %s" % workflow)
                return

            # On failure, parse errors and surface helpful diagnostics
            output_text = "".join(combined_output)
            failed_rules = self.parse_snakemake_errors(output_text)

            if not failed_rules:
                logging.error("Snakemake failed, but no rule-specific errors were parsed. See output above.")
                sys.exit(1)

            unique_logs = []
            seen = set()
            for rule in failed_rules:
                logs = self.find_logs_for_rule(rule, workflow, self.output, workflow_identifier)
                if logs:
                    for lp in logs:
                        if rule == "das_tool":
                            with open(lp) as f:
                                if "No bins were found, so DAS_tool cannot be run." in f.read():
                                    logging.info("--- Aviary -----------------------------------------------------------")
                                    logging.warning("No bins were found by any binners.")
                                    sys.exit(0)

                        logging.error(f"[{workflow_identifier}] Rule failed: {rule}; log: {lp}")
                        if lp not in seen:
                            unique_logs.append(lp)
                            seen.add(lp)
                else:
                    logs_root = self.logs_dir_root_for_workflow(self.output, workflow)
                    logging.error(f"[{workflow_identifier}] Rule failed: {rule}; no log files found under {logs_root}")

            # Dump log files to stderr
            for lp in unique_logs:
                try:
                    with open(lp, "r", encoding="utf-8", errors="replace") as fh:
                        sys.stderr.write(f"\n===== BEGIN LOG ({workflow_identifier}): {lp} =====\n")
                        sys.stderr.write(fh.read())
                        sys.stderr.write(f"\n===== END LOG ({workflow_identifier}): {lp} =====\n")
                except FileNotFoundError:
                    sys.stderr.write(f"\n===== LOG NOT FOUND ({workflow_identifier}): {lp} =====\n")
                except Exception as e:
                    sys.stderr.write(f"\n===== ERROR READING LOG ({workflow_identifier}) {lp}: {e} =====\n")
                sys.stderr.flush()

            sys.exit(1)

def fraction_to_percent(val):
    val = float(val)
    if val <= 1:
        return val * 100
    return val
