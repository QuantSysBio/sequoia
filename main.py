import os
import sys
import subprocess
import yaml
import shutil


def prepare_working_directory(server_project_dir_path, sequoia_project_dir_path):
    if not os.path.exists(sequoia_project_dir_path):
        os.mkdir(sequoia_project_dir_path)

    if not os.path.exists(f'{sequoia_project_dir_path}/data'):
        os.mkdir(f'{sequoia_project_dir_path}/data')
    
    if not os.path.exists(f'{sequoia_project_dir_path}/outfiles'):
        os.mkdir(f'{sequoia_project_dir_path}/outfiles')

    assert os.path.exists(f'{server_project_dir_path}/data/reads_fq')
    if not os.path.exists(f'{sequoia_project_dir_path}/data/reads_fq'):
        os.symlink(f'{server_project_dir_path}/data/reads_fq',
                f'{sequoia_project_dir_path}/data/reads_fq')

    shutil.copyfile('default_config.yaml',
                    sequoia_project_dir_path)


def combine_config_with_default(config):
    with open('default_config.yaml', 'r') as yamlfile:
        default_config = yaml.load(yamlfile, Loader=yaml.FullLoader)
    config = default_config.update(config)


def get_snakefile_path(config):
    return 'src/Snakefile'


if __name__ == '__main__':
    sequoia_dir_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    server_project_dir_path = sys.argv[1]

    project_name = os.path.basename(os.path.normpath(server_project_dir_path))
    sequoia_project_dir_path = f'projects/{project_name}'

    with open(f'{server_project_dir_path}/sequoia_config.yaml', 'r') as yamlfile:
        config = yaml.load(yamlfile, Loader=yaml.FullLoader)

    combine_config_with_default(config)
    prepare_working_directory(server_project_dir_path, sequoia_project_dir_path)
    snakefile_path = get_snakefile_path(config)
    slurm_config = config['slurm']

    p = subprocess.Popen([
        'snakemake', '--use-singularity', '--use-conda', 
        '--snakefile', snakefile_path,
        '--conda-prefix', f'{sequoia_dir_path}/.conda',
        '--cluster',
        f'sbatch -p {slurm_config["partition"]}',
        f'-N {slurm_config["nodes"]}',
        f'-c {slurm_config["ncpus"]}',
        f'--mem {slurm_config["mem"]}',
        f'-t {slurm_config["time"]}',
        f'--job-name {slurm_config["job_name"]}',
        f'-o {slurm_config["output"]}',
        f'-D {sequoia_project_dir_path}',
        '--cluster-cancel "scancel"',
        '--exclusive',
        '--conda-frontend', 'mamba',
        '-j', '3', #TODO
        '-w', '600',
        '--restart-times', '3',
        '--resources', 'load=100'],
        cwd=sequoia_project_dir_path)
    
    p.wait()
