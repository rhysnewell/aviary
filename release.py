#!/usr/bin/env python3

import io
from os.path import dirname, join
import extern

def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]

if __name__ == "__main__":
    version = get_version('aviary/__init__.py')
    print("version is {}".format(version))

    print("Updating dependencies ..")
    extern.run('pixi run --manifest-path aviary/pixi.toml update_deps')

    print("Checking if repo is clean ..")
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    extern.run('git tag v{}'.format(version))
    print("Now run 'git push && git push --tags' and GitHub actions will build and upload to PyPI")