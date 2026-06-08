#!/usr/bin/env python3
import os
from bird_tool_utils.pixi_deps import build_dep_defs

base_dir = os.path.join(os.path.dirname(__file__), "..")

build_dep_defs(
    pixi_toml_path=os.path.join(base_dir, "aviary", "pixi.toml"),
    pixi_lock_path=os.path.join(base_dir, "aviary", "pixi.lock"),
    environment_yml_path=os.path.join(base_dir, "admin", "environment.yml"),
    requirements_txt_path=os.path.join(base_dir, "admin", "requirements.txt"),
    conda_to_pip_name={"bird_tool_utils_python": "bird_tool_utils"},
    project_name="aviary-genome",  # conda name differs from workspace name "aviary"
)
