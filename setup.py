from os.path import dirname
from setuptools import setup, find_packages
import os


# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# We read requirements from the requirements.txt file, because that can be
# auto-generated from the pixi toml file.
base_dir = dirname(__file__)
with open(os.path.join(base_dir, "admin/requirements.txt")) as f:
    install_requires = f.read().splitlines()

# version= is deliberately absent: setuptools_scm derives it from the git tag
# (see [tool.setuptools_scm] in pyproject.toml). The former get_version() helper
# scraped __version__ out of aviary/__init__.py, which is now itself generated
# from the same source of truth -- keeping the scraper would have made the
# build read a file the build writes.
setup(
    name='aviary-genome',
    url='https://github.com/rhysnewell/aviary',
    license='GPL-3.0',
    author='Rhys Newell',
    author_email='rhys.newell94@gmail.com',
    description='aviary - metagenomics pipeline using long and short reads',
    long_description=long_description,
    long_description_content_type='text/markdown',
    zip_safe=False,
    packages=find_packages(),
    package_data={
        '': ['aviary/*'],
    },
    data_files=[(".", ["README.md", "LICENSE"])],
    include_package_data=True,
    install_requires= install_requires,
    entry_points={
          'console_scripts': [
              'aviary = aviary.aviary:main'
          ]
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
