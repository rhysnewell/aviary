import importlib.resources

def pixi_run_func():
    with importlib.resources.path("aviary", "pixi.toml") as manifest_path:
        return f"pixi run --frozen --manifest-path {manifest_path}"

pixi_run = pixi_run_func()
