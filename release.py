#!/usr/bin/env python3

import extern
import re
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the version number as an argument e.g. 0.13.0")
        sys.exit(1)
    version = sys.argv[1]

    yes_no = input(
        "Did you update CHANGELOG.md?\n\n"
    )
    if yes_no != "y":
        raise Exception("Please update the CHANGELOG.md file")

    yes_no = input(
        "Have all tests been run?\n\n"
    )
    if yes_no != "y":
        raise Exception("Please run all tests first")

    print("version is {}".format(version))

    print("Updating dependencies ..")
    extern.run('pixi run --manifest-path aviary/pixi.toml update_deps')

    # Update version in aviary/__init__.py
    print("Updating version in aviary/__init__.py")
    init_path = 'aviary/__init__.py'
    with open(init_path) as f:
        init_content = f.read()
    init_content = re.sub(
        r'^__version__ = ".*"',
        '__version__ = "{}"'.format(version),
        init_content,
        count=1,
        flags=re.MULTILINE,
    )
    with open(init_path, 'w') as f:
        f.write(init_content)

    print(
        "Checking for unexpected changes. If this fails you need to remove the git tag with 'git tag -d v{}'".format(version)
    )
    extern.run("if git diff --name-only | grep -qv 'aviary/__init__.py'; then echo 'Unexpected changed files:'; git diff --name-only; exit 1; fi")

    print("Committing the version file")
    extern.run('git commit -a -m "v{}"'.format(version))

    print("Tagging the release as v{}".format(version))
    extern.run('git tag v{}'.format(version))

    print("Now run 'git push && git push --tags' and GitHub actions will build and upload to PyPI")
