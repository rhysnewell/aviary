import json
import os
import sys
import signal
import subprocess

def handler(signum, frame):
     raise IOError
signal.signal(signal.SIGALRM, handler)


"""
Load the reference package. This will fail if the directory doesn't exist.
"""
def get_gtdb_path():
    try:
        GENERIC_PATH = os.environ['GTDBTK_DATA_PATH']
        return GENERIC_PATH
    except KeyError:
        print('\n' + '=' * 80)
        print(' ERROR '.center(80))
        print('_' * 80 + '\n')
        print("The 'GTDBTK_DATA_PATH' environment variable is not defined.".center(80) + '\n')
        print('Please set this variable to your reference data package.'.center(80))
        print('https://github.com/Ecogenomics/GTDBTk#installation'.center(80))
        print('Alternatively, use --gtdb-path flag.'.center(80))    
        print('=' * 80)

        signal.alarm(20)
        os.environ['GTDBTK_DATA_PATH'] = input('Input GTDBTK_DATA_PATH now:')
        try:
            subprocess.Popen('mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export GTDBTK_DATA_PATH=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset GTDBTK_DATA_PATH" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
                    (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['GTDBTK_DATA_PATH'], os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
        except KeyError:
            subprocess.Popen(
                'echo "export GTDBTK_DATA_PATH=%s" >> ~/.bashrc' %
                (os.environ['GTDBTK_DATA_PATH']), shell=True).wait()
        signal.alarm(0)


"""
Set the default enrichm db path.
"""
def get_enrichm_db_path():
    try:
        CONDA_PATH = os.environ['ENRICHM_DB']
        return CONDA_PATH
    except KeyError:
        print('\n' + '=' * 80)
        print(' ERROR '.center(80))
        print('_' * 80 + '\n')
        print("The 'ENRICHM_DB' environment variable is not defined.".center(80) + '\n')
        print('Please set this variable to your default server/home directory conda environment path.'.center(80))
        print('Alternatively, use --enrichm-db-path flag.'.center(80))
        print('=' * 80)
        signal.alarm(20)
        os.environ['ENRICHM_DB'] = input('Input ENRICHM_DB now:')
        try:
            subprocess.Popen(
                'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export ENRICHM_DB=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset ENRICHM_DB" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
                (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['ENRICHM_DB'],
                 os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
        except KeyError:
            subprocess.Popen(
                'echo "export ENRICHM_DB=%s" >> ~/.bashrc' %
                (os.environ['ENRICHM_DB']), shell=True).wait()
        signal.alarm(0)

"""
Set the default conda environment path.
"""
def get_conda_path():
    try:
        CONDA_PATH = os.environ['CONDA_ENV_PATH']
        return CONDA_PATH
    except KeyError:
        print('\n' + '=' * 80)
        print(' ERROR '.center(80))
        print('_' * 80 + '\n')
        print("The 'CONDA_ENV_PATH' environment variable is not defined.".center(80) + '\n')
        print('Please set this variable to your default server/home directory conda environment path.'.center(80))
        print('Alternatively, use --conda-prefix flag.'.center(80))    
        print('=' * 80)
        signal.alarm(20)
        os.environ['CONDA_ENV_PATH'] = input('Input CONDA_ENV_PATH now:')
        try:
            subprocess.Popen(
                'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export CONDA_ENV_PATH=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset CONDA_ENV_PATH" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
                (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['CONDA_ENV_PATH'],
                 os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
        except KeyError:
            subprocess.Popen(
                'echo "export CONDA_ENV_PATH=%s" >> ~/.bashrc' %
                (os.environ['CONDA_ENV_PATH']), shell=True).wait()
        signal.alarm(0)
