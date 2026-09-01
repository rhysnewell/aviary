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
from aviary.__init__ import (
    SHORT_READ_MAPPER_TO_COVERM, MAPPER_MODELS, MAPPERS_WITHOUT_MODELS,
    SHORT_READ_MODELS, LONG_READ_MODELS,
    DEFAULT_SHORT_READ_MAPPER, DEFAULT_LONG_READ_MAPPER,
    short_read_mapper_for_alignment,
)
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
from snakemake.common.configfile import load_configfile
from ruamel.yaml import YAML  # used for yaml reading with comments
from aviary import LONG_READ_TYPES, LONG_READ_TYPE_TO_SPADES, COVERAGE_JOB_STRATEGIES, COVERAGE_JOB_CUTOFF
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
# Subcommands that do not require short or long reads.
# Add new subcommand names here to allow them to run without reads.
SUBCOMMANDS_WITHOUT_READS = ['annotate', 'cluster', 'download_databases']

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
            self.skip_reads_check = args.skip_reads_check
        except AttributeError:
            self.skip_reads_check = any(w in SUBCOMMANDS_WITHOUT_READS for w in self.workflows)

        try:
            self.pplacer_threads = min(int(args.pplacer_threads), int(self.threads), 48)
        except AttributeError:
            self.pplacer_threads = min(int(self.threads), 48)

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
            self.semibin_mode = args.semibin_mode
            self.refinery_max_iterations = args.refinery_max_iterations
            self.refinery_max_retries = args.refinery_max_retries
            self.skip_abundances = args.skip_abundances
            self.skip_taxonomy = args.skip_taxonomy
            self.skip_singlem = args.skip_singlem
            self.filter_bins_min_completeness = args.min_completeness
            self.filter_bins_max_contamination = args.max_contamination
            if args.binning_only:
                self.skip_abundances = True
                self.skip_taxonomy = True
                self.skip_singlem = True
            self.binning_only = args.binning_only

            self.skip_binners = ["maxbin2", "concoct", "comebin", "taxvamb", "quickbin"]
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
                    elif binner == "quickbin":
                        self.skip_binners.remove("quickbin")
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
            self.semibin_mode = 'single'
            self.refinery_max_iterations = 5
            self.refinery_max_retries = 3
            self.skip_binners = ["none"]
            self.skip_abundances = False
            self.binning_only = False
            self.skip_taxonomy = False
            self.skip_singlem = False
            self.filter_bins_min_completeness = 50.0
            self.filter_bins_max_contamination = 5.0

        try:
            self.assembly = args.assembly
        except AttributeError:
            self.assembly = 'none'

        try:
            if args.host_filter != ['none']:
                self.host_filter = [os.path.abspath(ref_fil) for ref_fil in args.host_filter if ref_fil != 'none']
            else:
                self.host_filter = ['none']

            if args.gold_standard != ['none']:
                self.gold_standard = [os.path.abspath(p) for p in args.gold_standard]
            else:
                self.gold_standard = ['none']
            
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
            self.gold_standard = ['none']
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
            self.long_read_assembler = getattr(args, "long_read_assembler", "myloasm")
            self.long_read_mapper = getattr(args, "long_read_mapper", None)
            self.short_read_mapper = getattr(args, "short_read_mapper", None)
            self.long_read_mapper_model = getattr(args, "long_read_mapper_model", None)
            self.short_read_mapper_model = getattr(args, "short_read_mapper_model", None)
            self.minibwa_params = getattr(args, "minibwa_params", None)
            self.bwa_params = getattr(args, "bwa_params", None)
            self.strobealign_params = getattr(args, "strobealign_params", None)
            self.minimap2_params = getattr(args, "minimap2_params", None)
            self.rammap_params = getattr(args, "rammap_params", None)
        except AttributeError:
            self.longread_type = 'none'
            self.medaka_model = 'none'
            self.long_read_assembler = 'myloasm'
            self.long_read_mapper = None
            self.short_read_mapper = None
            self.long_read_mapper_model = None
            self.short_read_mapper_model = None
            self.minibwa_params = None
            self.bwa_params = None
            self.strobealign_params = None
            self.minimap2_params = None
            self.rammap_params = None
        self.guppy_model = getattr(args, 'guppy_model', 'r941_min_hac_g507')

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

        # Runs after the read inputs above are resolved, because the mapper
        # flags are validated against which reads were actually supplied.
        self._validate_mapper_selection()

        try:
            self.kmer_sizes = args.kmer_sizes
            self.use_unicycler = args.use_unicycler
            self.use_megahit = args.use_megahit
            self.coassemble = args.coassemble
            self.min_cov_long = args.min_cov_long
            self.min_cov_short = args.min_cov_short
            self.exclude_contig_cov = args.exclude_contig_cov
            self.exclude_contig_size = args.exclude_contig_size
            self.long_contig_size = args.include_contig_size
        except AttributeError:
            self.kmer_sizes = ['auto']
            self.use_unicycler = False
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
            min_completeness = args.min_completeness
            if min_completeness == 'none':
                self.min_completeness = ' '
            else:
                self.min_completeness = f'--min-completeness {min_completeness}'

            max_contamination = args.max_contamination
            if max_contamination == 'none':
                self.max_contamination = ' '
            else:
                self.max_contamination = f'--max-contamination {max_contamination}'

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


    def _validate_mapper_selection(self):
        """Reject impossible --*-mapper / --*-mapper-model / --minibwa-params
        combinations before the run starts.

        These all used to surface only once a mapping job actually ran, i.e.
        after assembly had already burned hours of walltime -- or, worse, not at
        all: CoverM accepts a short-read preset for long reads and returns a
        well-formed table of near-zero depths rather than failing, so the
        mistake reads as bad numbers instead of an error.

        Requires self.pe1/self.pe2/self.longreads to be resolved already.
        """
        has_short_reads = self.pe1 != 'none' and self.pe1
        has_long_reads = self.longreads != 'none' and self.longreads

        # A mapper is only meaningful alongside the reads it would map. These
        # are None unless the user passed the flag (see aviary.py), so an
        # unrelated default never trips this.
        for mapper_value, model_value, reads_present, mapper_flag, model_flag, reads_flag in (
            (self.short_read_mapper, self.short_read_mapper_model, has_short_reads,
             "--short-read-mapper", "--short-read-mapper-model", "-1/-2/--interleaved/--coupled"),
            (self.long_read_mapper, self.long_read_mapper_model, has_long_reads,
             "--long-read-mapper", "--long-read-mapper-model", "-l/--longreads"),
        ):
            if reads_present:
                continue
            for value, flag in ((mapper_value, mapper_flag), (model_value, model_flag)):
                if value is not None:
                    logging.error(
                        f"{flag} was given as {value!r}, but no reads were supplied for it "
                        f"to map. Provide reads with {reads_flag}, or drop {flag}."
                    )
                    sys.exit(-1)

        # Resolve to the documented defaults now that "was it given?" has been
        # answered. Everything downstream sees a concrete mapper name.
        if self.short_read_mapper is None:
            self.short_read_mapper = DEFAULT_SHORT_READ_MAPPER
        if self.long_read_mapper is None:
            self.long_read_mapper = DEFAULT_LONG_READ_MAPPER

        # strobealign-aemb shells out to `strobealign --aemb` directly and
        # never produces a BAM file (see get_coverage.py); it only ever
        # yields a coverage table. concoct/semibin/comebin/quickbin all read
        # data/binning_bams/ directly and cannot function without one,
        # regardless of semibin-mode -- they would otherwise fail deep in the
        # DAG with an opaque "Expected file '*.bam' does not exist", with no
        # indication the root cause is the mapper choice. rosella/metabat2/
        # vamb/taxvamb/maxbin2 read the coverage table instead and are
        # unaffected. Auto-skip the BAM-only binners with a warning rather
        # than fail mid-run.
        if self.short_read_mapper == "strobealign-aemb":
            bam_only_binners = ("semibin", "concoct", "comebin", "quickbin")
            newly_skipped = [b for b in bam_only_binners if b not in self.skip_binners]
            if newly_skipped:
                logging.warning(
                    f"--short-read-mapper strobealign-aemb produces no BAM file, so the "
                    f"following binners cannot run without one. Automatically skipping: "
                    f"{', '.join(newly_skipped)}."
                )
                self.skip_binners.extend(newly_skipped)

        # --*-mapper-model is only meaningful for mappers with more than one
        # CoverM preset (minimap2, rammap). Validated here rather than via
        # argparse choices=, since the valid set depends on the paired mapper
        # flag -- and on read length, because CoverM exposes both short- and
        # long-read presets through the same -p option.
        for model_value, mapper_value, flag_name, valid_for_length, length_name in (
            (self.short_read_mapper_model, self.short_read_mapper,
             "--short-read-mapper-model", SHORT_READ_MODELS, "short"),
            (self.long_read_mapper_model, self.long_read_mapper,
             "--long-read-mapper-model", LONG_READ_MODELS, "long"),
        ):
            if model_value is None:
                continue
            if mapper_value in MAPPERS_WITHOUT_MODELS:
                logging.error(
                    f"{flag_name} was given but {mapper_value!r} has no selectable model."
                )
                sys.exit(-1)
            if mapper_value in MAPPER_MODELS and model_value not in MAPPER_MODELS[mapper_value]:
                logging.error(
                    f"Wrong model chosen for mapper {mapper_value!r}: {model_value!r} is not valid. "
                    f"Valid models: {MAPPER_MODELS[mapper_value]}"
                )
                sys.exit(-1)
            if model_value not in valid_for_length:
                logging.error(
                    f"{flag_name} was given as {model_value!r}, which is not a {length_name}-read "
                    f"preset. CoverM would accept it and return a well-formed table of near-zero "
                    f"depths rather than failing. Valid {length_name}-read models: "
                    f"{list(valid_for_length)}"
                )
                sys.exit(-1)

        if self.minibwa_params is not None and "minibwa" not in (self.short_read_mapper, self.long_read_mapper):
            logging.error(
                "--minibwa-params was given but minibwa is not selected as "
                "--short-read-mapper or --long-read-mapper."
            )
            sys.exit(-1)

        # Raw CoverM per-aligner passthroughs, same "given without the mapper
        # it applies to" guard as --minibwa-params above. bwa-mem/bwa-mem2 and
        # strobealign are short-read-only mappers, so --bwa-params/
        # --strobealign-params are only checked against --short-read-mapper;
        # minimap2/rammap are valid for either read length, so --minimap2-params/
        # --rammap-params are checked against both. strobealign-aemb
        # deliberately does not count for --strobealign-params: it shells out
        # to `strobealign --aemb` directly inside CoverM and does not accept it.
        if self.bwa_params is not None and self.short_read_mapper not in ("bwa-mem", "bwa-mem2"):
            logging.error(
                "--bwa-params was given but --short-read-mapper is not bwa-mem or bwa-mem2."
            )
            sys.exit(-1)
        if self.strobealign_params is not None and self.short_read_mapper != "strobealign":
            logging.error(
                "--strobealign-params was given but --short-read-mapper is not strobealign "
                "(strobealign-aemb does not accept it)."
            )
            sys.exit(-1)
        if self.minimap2_params is not None and "minimap2" not in (self.short_read_mapper, self.long_read_mapper):
            logging.error(
                "--minimap2-params was given but minimap2 is not selected as "
                "--short-read-mapper or --long-read-mapper."
            )
            sys.exit(-1)
        if self.rammap_params is not None and "rammap" not in (self.short_read_mapper, self.long_read_mapper):
            logging.error(
                "--rammap-params was given but rammap is not selected as "
                "--short-read-mapper or --long-read-mapper."
            )
            sys.exit(-1)


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
            if len(self.assembly) > 1 and self.semibin_mode != "multi":
                logging.error(
                    "Multiple assemblies provided but --semibin-mode is not 'multi'. "
                    "Pass --semibin-mode multi to enable SemiBin2 multi-sample binning, "
                    "or provide a single pre-concatenated assembly."
                )
                sys.exit(-1)
            if len(self.assembly) == 1 and self.semibin_mode == "multi":
                logging.warning(
                    "--semibin-mode multi was requested but only one assembly was provided. "
                    "SemiBin2 will run but multi-sample binning requires at least two assemblies "
                    "to be meaningful. Pass multiple assemblies with --assembly A.fasta B.fasta "
                    "to use multi-sample mode properly."
                )
        elif self.assembly is None:
            if self.semibin_mode == "multi":
                logging.error(
                    "--semibin-mode multi requires at least two assemblies passed with "
                    "--assembly A.fasta B.fasta, but no assembly was provided. Reads-only "
                    "runs assemble a single set of contigs, which multi-sample binning "
                    "cannot use. Provide multiple assemblies, or drop --semibin-mode multi."
                )
                sys.exit(-1)
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
        conf["skip_reads_check"] = self.skip_reads_check
        conf["gsa"] = self.gold_standard
        conf["gsa_mappings"] = self.gsa_mappings
        conf["skip_binners"] = self.skip_binners
        conf["skip_abundances"] = self.skip_abundances
        conf["skip_taxonomy"] = self.skip_taxonomy
        conf["skip_singlem"] = self.skip_singlem
        conf["binning_only"] = self.binning_only
        conf["semibin_model"] = self.semibin_model
        conf["semibin_mode"] = self.semibin_mode
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
        # spades_assembly.py accepts a narrower --long-read-type vocabulary
        # than aviary's own (see LONG_READ_TYPE_TO_SPADES) -- resolved here,
        # not passed raw, so rs/sq/ccs/hifi don't crash it with "invalid choice".
        # .get(..., "none") rather than a bare lookup: self.longread_type is
        # the literal string "none" whenever no long reads were supplied at
        # all (annotate, short-read-only recover, ...), and
        # LONG_READ_TYPE_TO_SPADES only maps the six real read types.
        # spades_assembly.py is only ever invoked on the hybrid/long-read
        # assembly path, which is unreachable when there are no long reads,
        # so this fallback value is never actually read -- it only needs to
        # exist because the .smk rules interpolate config values
        # unconditionally, same reason every other inapplicable field here
        # (long_read_mapper_model, minibwa_params, ...) is written as "none"
        # rather than omitted.
        conf["long_read_type_spades"] = LONG_READ_TYPE_TO_SPADES.get(self.longread_type, "none")
        conf["long_read_assembler"] = self.long_read_assembler
        conf["long_read_mapper"] = self.long_read_mapper
        conf["long_read_mapper_model"] = self.long_read_mapper_model or "none"
        conf["minibwa_params"] = self.minibwa_params or "none"
        conf["bwa_params"] = self.bwa_params or "none"
        conf["strobealign_params"] = self.strobealign_params or "none"
        conf["minimap2_params"] = self.minimap2_params or "none"
        conf["rammap_params"] = self.rammap_params or "none"
        # Resolved to the value CoverM's -p expects, so the Snakefiles and the
        # coverage scripts can use it verbatim. Long reads stay as the family
        # name since their preset depends on --long-read-type (or
        # --long-read-mapper-model, resolved downstream in the scripts).
        if self.short_read_mapper in MAPPER_MODELS and self.short_read_mapper_model is not None:
            conf["short_read_mapper"] = f"{self.short_read_mapper}-{self.short_read_mapper_model}"
        else:
            conf["short_read_mapper"] = SHORT_READ_MAPPER_TO_COVERM[self.short_read_mapper]
        # `coverm genome` (per-genome relative abundance) cannot run
        # strobealign-aemb -- it is `coverm contig`-only. get_abundances.py and
        # the raw dereplicate_and_get_abundances_* rules use this key instead
        # of short_read_mapper so they always get a real -p value.
        conf["short_read_mapper_aligner"] = short_read_mapper_for_alignment(conf["short_read_mapper"])
        conf["medaka_model"] = self.medaka_model
        conf["guppy_model"] = self.guppy_model
        conf["kmer_sizes"] = self.kmer_sizes
        conf["use_unicycler"] = self.use_unicycler
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
        # conf["strain_analysis"] = self.strain_analysis
        conf["checkm2_db_folder"] = self.checkm2_db
        conf["use_checkm2_scores"] = self.use_checkm2_scores
        conf["mag_directory"] = self.mag_directory
        conf["mag_extension"] = self.mag_extension
        conf["previous_runs"] = self.previous_runs
        conf["min_completeness"] = self.min_completeness
        conf["max_contamination"] = self.max_contamination
        conf["filter_bins_min_completeness"] = self.filter_bins_min_completeness
        conf["filter_bins_max_contamination"] = self.filter_bins_max_contamination
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

        if dryrun and self.longreads != 'none':
            if self.long_read_assembler == 'flye':
                print("logs/flye_assembly.log")
            elif self.long_read_assembler == 'myloasm':
                print("logs/myloasm_assembly.log")
            else:
                raise Exception("Programming error: unexpected long_read_assembler value.")

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
                continue

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
