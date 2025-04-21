import io
from os.path import dirname, join
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

def get_version(relpath):
  """Read version info from a file without importing it"""
  for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
    if "__version__" in line:
      if '"' in line:
        # __version__ = "0.9"
        return line.split('"')[1]
      elif "'" in line:
        return line.split("'")[1]


setup(
    name='aviary-genome',
    version=get_version("aviary/__init__.py"),
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
