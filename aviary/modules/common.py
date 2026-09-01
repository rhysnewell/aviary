import importlib.resources
import os
import time

def pixi_run_func():
    with importlib.resources.path("aviary", "pixi.toml") as manifest_path:
        # Do not pub --frozen here, otherwise dev becomes more confusing because
        # pixi.toml changes have no effect.
        return f"pixi run --manifest-path {manifest_path}"

pixi_run = pixi_run_func()

# A unique identifier for the current workflow invocation. This is used
# when creating log directories so that retries from the same workflow do
# not overwrite previous log files.
workflow_identifier = time.strftime("%Y%m%d_%H%M%S")

def setup_log(log_dir_base: str, attempt: int) -> str:
    """Return a unique log file path for a given rule attempt.
    Parameters
    ----------
    log_dir_base: str
        Directory in which log files for a rule should be stored.  This
        function will create a subdirectory named with the workflow
        identifier and place attempt specific log files inside it.
    attempt: int
        The Snakemake retry attempt number.
    Returns
    -------
    str
        Path to a log file unique to the workflow invocation and attempt.
    """

    log_dir = os.path.join(log_dir_base, workflow_identifier)
    os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, f"attempt{attempt}.log")

def get_semibin_mode(config: dict) -> str:
    """Return the configured SemiBin mode.
    Supports the older semibin_multi config key for compatibility with early test configs.
    """
    if "semibin_mode" in config:
        return config["semibin_mode"]
    return "multi" if config.get("semibin_multi", False) else "single"

def primary_fasta(config: dict) -> str:
    """A single representative assembly path for whole-run QC/reporting.

    config["fasta"] is a list of paths when SemiBin2 multi-sample binning is
    active (config["semibin_mode"] == "multi", one entry per --assembly
    file), and a single path string otherwise. QC/reporting rules need one
    file, not the raw list. In multi mode, the SemiBin2-concatenated fasta
    (the file binning itself works from) is the correct single reference --
    not any one of the pre-concatenation inputs -- since it's the assembly
    actually used downstream.
    """
    if get_semibin_mode(config) == "multi":
        return "data/semibin_multi_prep/concatenated.fa"
    return config["fasta"]
