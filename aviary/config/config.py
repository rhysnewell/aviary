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
    try:
        with open(format('%s/etc/conda/activate.d/aviary.sh' % os.environ['CONDA_PREFIX'])) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                # if 'export' not in line:
                #     continue
                # Remove leading `export `, if you have those
                # then, split name / value pair
                # key, value = line.replace('export ', '', 1).strip().split('=', 1)
                try:
                    key, value = line.strip().split('=', 1)
                    key = key.strip('export ')
                    # print(key, value)
                    os.environ[key] = value  # Load to local environ
                except ValueError:
                    continue
    except FileNotFoundError:
        # File not found so going to have to create it
        pass

def source_bashrc():
    try:
        with open('%s/.bashrc' % os.environ['HOME']) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                # if 'export' not in line:
                #     continue
                # Remove leading `export `, if you have those
                # then, split name / value pair
                # key, value = line.replace('export ', '', 1).strip().split('=', 1)
                try:
                    key, value = line.strip().split('=', 1)
                    key = key.strip('export ')
                    os.environ[key] = value  # Load to local environ
                except ValueError:
                    continue
    except FileNotFoundError:
        # File not found so going to have to create it
        pass

"""
Load the reference package. This will fail if the directory doesn't exist.
"""
def get_software_db_path(db_name='CONDA_ENV_PATH', software_flag='--conda-prefix'):
    try:
        source_conda_env()
        SW_PATH = os.environ[db_name]
        return SW_PATH
    except KeyError:
        try:
            source_conda_env()
            CONDA_PATH = os.environ[db_name]
            return CONDA_PATH
        except KeyError:
            try:
                source_bashrc()
                CONDA_PATH = os.environ[db_name]
                return CONDA_PATH
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
                try:
                    conda_prefix = os.environ['CONDA_PREFIX']
                    # make the directory
                    os.makedirs(f"{conda_prefix}/etc/conda/activate.d/", exist_ok=True)
                    os.makedirs(f"{conda_prefix}/etc/conda/deactivate.d/", exist_ok=True)
                    # add the export to the activate script
                    with open(f"{conda_prefix}/etc/conda/activate.d/aviary.sh", 'a') as f:
                        f.write(f'export {db_name}={os.environ[db_name]}\n')
                    
                    # add the unset to the deactivate script
                    with open(f"{conda_prefix}/etc/conda/deactivate.d/aviary.sh", 'a') as f:
                        f.write(f'unset {db_name}\n')

                except KeyError:
                    # put the export in the bashrc
                    with open(f"{os.environ['HOME']}/.bashrc", 'a') as f:
                        f.write(f'export {db_name}={os.environ[db_name]}\n')

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
    try:
        conda_prefix = os.environ['CONDA_PREFIX']
        # make the directory
        os.makedirs(f"{conda_prefix}/etc/conda/activate.d/", exist_ok=True)
        os.makedirs(f"{conda_prefix}/etc/conda/deactivate.d/", exist_ok=True)
        # add the export to the activate script
        with open(f"{conda_prefix}/etc/conda/activate.d/aviary.sh", 'a') as f:
            f.write(f'export {db_name}={os.environ[db_name]}\n')
        
        # add the unset to the deactivate script
        with open(f"{conda_prefix}/etc/conda/deactivate.d/aviary.sh", 'a') as f:
            f.write(f'unset {db_name}\n')
    except KeyError:
        # put the export in the bashrc
        with open(f"{os.environ['HOME']}/.bashrc", 'a') as f:
            f.write(f'export {db_name}={os.environ[db_name]}\n')
