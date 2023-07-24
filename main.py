import os
import sys
import subprocess
import yaml
import shutil


def prepare_working_directory(server_project_dir_path, sequoia_project_dir_path, sequoia_dir_path, config):
    # main directory
    if not os.path.exists(sequoia_project_dir_path):
        os.makedirs(sequoia_project_dir_path)

    # reads
    if not os.path.exists(f'{sequoia_project_dir_path}/data'):
        os.mkdir(f'{sequoia_project_dir_path}/data')

    assert os.path.exists(f'{server_project_dir_path}/data/reads_fq')
    if not os.path.exists(f'{sequoia_project_dir_path}/data/reads_fq'):
        shutil.copytree(f'{server_project_dir_path}/data/reads_fq',
                        f'{sequoia_project_dir_path}/data/reads_fq')
    
    # outfiles
    if not os.path.exists(f'{sequoia_project_dir_path}/outfiles'):
        os.mkdir(f'{sequoia_project_dir_path}/outfiles')

    # references
    if not os.path.exists(f'{sequoia_project_dir_path}/references'):
        shutil.copytree(f'{sequoia_dir_path}/references',
                        f'{sequoia_project_dir_path}/references')
        
    #config
    with open(f'{sequoia_project_dir_path}/config.yaml', 'w') as yamlfile:
        yaml.dump(config, yamlfile)


def _combine_config_dict(config, default_config):
    for key, item in default_config.items():
        if key not in config:
            config[key] = item
        elif isinstance(item, dict):
            config[key] = _combine_config_dict(config[key], default_config[key])
    
    return config


def combine_config_with_default(config):
    with open('default_config.yaml', 'r') as yamlfile:
        default_config = yaml.load(yamlfile, Loader=yaml.FullLoader)
    return _combine_config_dict(config, default_config)


def get_snakefile_path():
    return os.path.abspath('src/Snakefile')


def return_results(server_project_dir_path, sequoia_project_dir_path):
    sequoia_results_dir_path = f'{sequoia_project_dir_path}/results'
    server_results_dir_path = f'{server_project_dir_path}/sequoia_results'
    
    os.mkdir(server_results_dir_path)    
    shutil.move(f'{sequoia_results_dir_path}/tr_2_prot/proteome_expressed.fasta',
                f'{server_results_dir_path}/proteome_expressed.fasta')
    
    shutil.move(f'{sequoia_results_dir_path}/tr_2_prot/gene_tr_prot_SP.csv',
                f'{server_results_dir_path}/gene_tr_prot_SP.csv')

    for file_name in ['gene_expression_tpm.csv', 'transcript_expression_tpm.csv',
                      'gene_expression.csv', 'transcript_expression.csv']:
        shutil.move(f'{sequoia_results_dir_path}/transcriptome_quant/{file_name}',
                    f'{server_results_dir_path}/{file_name}')

    os.mkdir(f'{server_results_dir_path}/reports')
    shutil.move(f'{sequoia_results_dir_path}/multiqc/multiqc_report.html',
                f'{server_results_dir_path}/reports/multiqc_report.html')
    
    shutil.move(f'{sequoia_results_dir_path}/benchmarks',
                f'{server_results_dir_path}/reports/benchmarks')
    
    shutil.move(f'{sequoia_results_dir_path}/logs',
                f'{server_results_dir_path}/reports/logs')


if __name__ == '__main__':
    sequoia_dir_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    server_project_dir_path = sys.argv[1]

    project_name = os.path.basename(os.path.normpath(server_project_dir_path))
    sequoia_project_dir_path = os.path.abspath(f'projects/{project_name}')

    with open(f'{server_project_dir_path}/sequoia_config.yaml', 'r') as yamlfile:
        config = yaml.load(yamlfile, Loader=yaml.FullLoader)

    config = combine_config_with_default(config)
    prepare_working_directory(server_project_dir_path, sequoia_project_dir_path, sequoia_dir_path, config)
    snakefile_path = get_snakefile_path()
    slurm_config = config['slurm']

    p = subprocess.Popen([
        'snakemake', '--use-singularity', '--use-conda',
        '--snakefile', snakefile_path,
        '--cluster', '/data/software/slurm/bin/sbatch '
        f'-p {slurm_config["partition"]} '
        f'-N {slurm_config["nodes"]} '
        f'-c {slurm_config["ncpus"]} '
        f'--mem {slurm_config["mem"]} '
        f'-t {slurm_config["time"]} '
        f'--job-name {slurm_config["job_name"]} '
        f'-o {slurm_config["output"]} '
        f'-D {sequoia_project_dir_path} --exclusive',
        '--cluster-cancel', 'scancel',
        '--conda-frontend', 'mamba',
        '--singularity-args', f'-B {sequoia_dir_path}:{sequoia_dir_path}',
        '--conda-prefix', f'{sequoia_dir_path}/.snakemake',
        '--singularity-prefix', f'{sequoia_dir_path}/.snakemake',
        '-j', '3',
        '-w', '600',
        '--restart-times', '1',
        '--resources', 'load=100'],
        cwd=sequoia_project_dir_path)
    p.wait()

    return_results(server_project_dir_path, sequoia_project_dir_path)
    shutil.rmtree(sequoia_project_dir_path)
