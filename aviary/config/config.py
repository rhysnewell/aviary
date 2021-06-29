import os
import signal
import subprocess

"""
Function to handle signal IOErrors after missing input
"""
def handler(signum, frame):
     raise IOError
signal.signal(signal.SIGALRM, handler)


def source_conda_env():
    with open(format('%s/etc/conda/activate.d/aviary.sh' % os.environ['CONDA_PREFIX'])) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            # if 'export' not in line:
            #     continue
            # Remove leading `export `, if you have those
            # then, split name / value pair
            # key, value = line.replace('export ', '', 1).strip().split('=', 1)
            key, value = line.strip().split('=', 1)
            os.environ[key] = value  # Load to local environ

def source_bashrc():
    with open('%s/.bashrc' % os.environ['HOME']) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            # if 'export' not in line:
            #     continue
            # Remove leading `export `, if you have those
            # then, split name / value pair
            # key, value = line.replace('export ', '', 1).strip().split('=', 1)
            key, value = line.strip().split('=', 1)
            os.environ[key] = value  # Load to local environ

"""
Load the reference package. This will fail if the directory doesn't exist.
"""
def get_gtdb_path():
    try:
        source_conda_env()
        GTDB_PATH = os.environ['GTDBTK_DATA_PATH']
        return GTDB_PATH
    except KeyError:
        try:
            source_conda_env()
            GTDB_PATH = os.environ['GTDBTK_DATA_PATH']
            return GTDB_PATH
        except KeyError:
            try:
                source_bashrc()
                GTDB_PATH = os.environ['GTDBTK_DATA_PATH']
                return GTDB_PATH
            except KeyError:
                print('\n' + '=' * 100)
                print(' ERROR '.center(100))
                print('_' * 100 + '\n')
                print("The 'GTDBTK_DATA_PATH' environment variable is not defined.".center(100) + '\n')
                print('Please set this variable to your reference data package.'.center(100))
                print('https://github.com/Ecogenomics/GTDBTk#installation'.center(100))
                print('Alternatively, use --gtdb-path flag.'.center(100))
                print('=' * 100)

                signal.alarm(120)
                os.environ['GTDBTK_DATA_PATH'] = input('Input GTDBTK_DATA_PATH now:')
                try:
                    subprocess.Popen('mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export GTDBTK_DATA_PATH=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset GTDBTK_DATA_PATH" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
                            (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['GTDBTK_DATA_PATH'], os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
                except KeyError:
                    subprocess.Popen(
                        'echo "export GTDBTK_DATA_PATH=%s" >> ~/.bashrc' %
                        (os.environ['GTDBTK_DATA_PATH']), shell=True).wait()
                signal.alarm(0)
                print('=' * 100)
                print('Reactivate your aviary conda environment or source ~/.bashrc to suppress this message.'.center(100))
                print('=' * 100)

                return os.environ['GTDBTK_DATA_PATH']

def set_gtdb_path(path):
    os.environ['GTDBTK_DATA_PATH'] = path
    try:
        subprocess.Popen(
            'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export GTDBTK_DATA_PATH=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset GTDBTK_DATA_PATH" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
            (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['GTDBTK_DATA_PATH'],
             os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
    except KeyError:
        subprocess.Popen(
            'echo "export GTDBTK_DATA_PATH=%s" >> ~/.bashrc' %
            (os.environ['GTDBTK_DATA_PATH']), shell=True).wait()

"""
Set the default enrichm db path.
"""

def get_enrichm_db_path():
    try:
        source_conda_env()
        ENRICHM_DB = os.environ['ENRICHM_DB']
        return ENRICHM_DB
    except KeyError:
        try:
            source_conda_env()
            ENRICHM_DB = os.environ['ENRICHM_DB']
            return ENRICHM_DB
        except KeyError:
            try:
                source_bashrc()
                ENRICHM_DB = os.environ['ENRICHM_DB']
                return ENRICHM_DB
            except KeyError:
                print('\n' + '=' * 100)
                print(' ERROR '.center(100))
                print('_' * 100 + '\n')
                print("The 'ENRICHM_DB' environment variable is not defined.".center(100) + '\n')
                print('Please set this variable to your default server/home directory conda environment path.'.center(100))
                print('Alternatively, use --enrichm-db-path flag.'.center(100))
                print('=' * 100)
                signal.alarm(120)
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
                print('=' * 100)
                print('Reactivate your aviary conda environment or source ~/.bashrc to suppress this message.'.center(100))
                print('=' * 100)

                return os.environ['ENRICHM_DB']


def set_enrichm_db_path(path):
    os.environ['ENRICHM_DB'] = path
    try:
        subprocess.Popen(
            'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export ENRICHM_DB=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset ENRICHM_DB" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
            (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['ENRICHM_DB'],
             os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
    except KeyError:
        subprocess.Popen(
            'echo "export ENRICHM_DB=%s" >> ~/.bashrc' %
            (os.environ['ENRICHM_DB']), shell=True).wait()

"""
Set the default enrichm db path.
"""


def get_busco_db_path():
    try:
        source_conda_env()
        BUSCO_DB = os.environ['BUSCO_DB']
        return BUSCO_DB
    except KeyError:
        try:
            source_conda_env()
            BUSCO_DB = os.environ['BUSCO_DB']
            return BUSCO_DB
        except KeyError:
            try:
                source_bashrc()
                BUSCO_DB = os.environ['BUSCO_DB']
                return BUSCO_DB
            except KeyError:
                print('\n' + '=' * 100)
                print(' ERROR '.center(100))
                print('_' * 100 + '\n')
                print("The 'BUSCO_DB' environment variable is not defined.".center(100) + '\n')
                print('Please set this variable to your default server/home directory conda environment path now.'.center(100))
                print('Alternatively, use --busco-db-path flag.'.center(100))
                print('=' * 100)
                signal.alarm(120)
                os.environ['BUSCO_DB'] = input('Input BUSCO_DB now:')
                try:
                    subprocess.Popen(
                        'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export BUSCO_DB=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset BUSCO_DB" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
                        (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['BUSCO_DB'],
                         os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
                except KeyError:
                    subprocess.Popen(
                        'echo "export BUSCO_DB=%s" >> ~/.bashrc' %
                        (os.environ['BUSCO_DB']), shell=True).wait()
                signal.alarm(0)
                print('=' * 100)
                print('Reactivate your aviary conda environment or source ~/.bashrc to suppress this message.'.center(100))
                print('=' * 100)

                return os.environ['BUSCO_DB']


def set_busco_db_path(path):
    os.environ['BUSCO_DB'] = path
    try:
        subprocess.Popen(
            'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export BUSCO_DB=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset BUSCO_DB" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
            (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['BUSCO_DB'],
             os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
    except KeyError:
        subprocess.Popen(
            'echo "export BUSCO_DB=%s" >> ~/.bashrc' %
            (os.environ['BUSCO_DB']), shell=True).wait()


"""
Set the default conda environment path.
"""

def get_conda_path():
    try:
        source_conda_env()
        CONDA_PATH = os.environ['CONDA_ENV_PATH']
        return CONDA_PATH
    except KeyError:
        try:
            source_conda_env()
            CONDA_PATH = os.environ['CONDA_ENV_PATH']
            return CONDA_PATH
        except KeyError:
            try:
                source_bashrc()
                CONDA_PATH = os.environ['CONDA_ENV_PATH']
                return CONDA_PATH
            except KeyError:
                print('\n' + '=' * 100)
                print(' ERROR '.center(100))
                print('_' * 100 + '\n')
                print("The 'CONDA_ENV_PATH' environment variable is not defined.".center(100) + '\n')
                print('Please set this variable to your default server/home directory conda environment path.'.center(100))
                print('Alternatively, use --conda-prefix flag.'.center(100))
                print('=' * 100)
                signal.alarm(120)
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
                print('=' * 100)
                print('Reactivate your aviary conda environment or source ~/.bashrc to suppress this message.'.center(100))
                print('=' * 100)

                return os.environ['CONDA_ENV_PATH']


def set_conda_path(path):
    os.environ['CONDA_ENV_PATH'] = path
    try:
        subprocess.Popen(
            'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export CONDA_ENV_PATH=%s" >> %s/etc/conda/activate.d/aviary.sh; echo "unset CONDA_ENV_PATH" >> %s/etc/conda/deactivate.d/aviary.sh; ' %
            (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], os.environ['CONDA_ENV_PATH'],
             os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX']), shell=True).wait()
    except KeyError:
        subprocess.Popen(
            'echo "export CONDA_ENV_PATH=%s" >> ~/.bashrc' %
            (os.environ['CONDA_ENV_PATH']), shell=True).wait()