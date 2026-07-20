import os
import signal
import subprocess
from aviary.modules.common import pixi_run

"""
Function to handle signal IOErrors after missing input
"""
def handler(signum, frame):
     raise IOError
signal.signal(signal.SIGALRM, handler)


def configure_variable(variable, value):
    os.environ[variable] = value
    subprocess.run(f"{pixi_run} conda env config vars set {variable}={value}".split(), check=True, capture_output=True)
    try:
        subprocess.run(f"pixi run --frozen conda env config vars set {variable}={value}".split(), check=True, capture_output=True)
    except subprocess.CalledProcessError:
        subprocess.run(f"conda env config vars set {variable}={value}".split(), check=True, capture_output=True)

"""
Load the reference package. This will fail if the directory doesn't exist.
"""
def get_software_db_path(db_name='CONDA_ENV_PATH', software_flag='--conda-prefix'):
    try:
        SW_PATH = os.environ[db_name]
        return SW_PATH
    except KeyError:
        try:
            SW_PATH = os.environ[db_name]
            configure_variable(db_name, SW_PATH)
            return SW_PATH
        except KeyError:
            print('\n' + '=' * 100)
            print(' ERROR '.center(100))
            print('_' * 100 + '\n')
            print(f"The '{db_name}' environment variable is not defined.".center(100) + '\n')
            print(f'Please set this variable to your default server/home directory containing {db_name}.'.center(100))
            print(f'Alternatively, use {software_flag} flag.'.center(100))
            print(f'Note: This variable must point to the DIRECTORY containing the files, not the files themselves'.center(100))
            print('Note: This can be set to an arbitrary string if you do not need this database'.center(100))
            print('=' * 100)
            signal.alarm(120)
            os.environ[db_name] = input(f'Input path to directory for {db_name} now:').strip()
            configure_variable(db_name, os.environ[db_name])

            signal.alarm(0)
            print('=' * 100)
            # print('Reactivate your aviary conda environment or source ~/.bashrc to suppress this message.'.center(100))
            # print('=' * 100)

            return os.environ[db_name]


"""
Sets an environmental variable and appends it to the conda activation script
"""
def set_db_path(path, db_name='CONDA_ENV_PATH'):
    os.environ[db_name] = path.strip()
    configure_variable(db_name, path.strip())
