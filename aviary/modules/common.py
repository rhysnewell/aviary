import importlib.resources

def pixi_run_func():
    pixi_path = importlib.resources.path("aviary", "pixi.toml")
    return f"pixi run --frozen --manifest-path {pixi_path}"

pixi_run = pixi_run_func()
