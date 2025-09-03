import importlib.resources
import os
import time

def pixi_run_func():
    with importlib.resources.path("aviary", "pixi.toml") as manifest_path:
        return f"pixi run --frozen --manifest-path {manifest_path}"

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
