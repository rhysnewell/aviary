import os

def pixi_run_func():
    aviary_basedir = os.path.join(os.path.dirname(__file__), "..", "..")
    return f"pixi run --frozen --manifest-path {aviary_basedir}/pixi.toml"

pixi_run = pixi_run_func()
