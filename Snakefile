shell.prefix("set -euo pipefail;")
shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

container: "docker://snakemake/snakemake:v6.8.0"
min_version("6.0")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
if config["reference"]['ERCC_normalization']:
    include: snakefiles + 'pre_process_workflow.py'
if config['workflow']['isoform_assembly'] or config['workflow']['reference_quant']:
    include: snakefiles + "folders.py"
    include: snakefiles + "reads_wildcards.py"
    include: snakefiles + "pre_process.py"
if config['workflow']['isoform_assembly']:
    include: snakefiles + "STAR.py"
    include: snakefiles + "filter_BAM.py"
    include: snakefiles + "stringtie.py"
    include: snakefiles + "salmon_extended.py"
    include: snakefiles + "transcriptome_2_proteome_extended.py"
    include: snakefiles + "multiqc_extended.py"
if config['workflow']['reference_quant']:
    include: snakefiles + 'salmon_refrence.py'
    include: snakefiles + 'transcriptome_2_proteome_refrence.py'
    include: snakefiles + 'multiqc_refrence.py'
if config['workflow']['exhaustive_ORFs']:
    include: snakefiles + 'exhaustive_workflow.py'

qc_path = []

if config['workflow']['isoform_assembly']:
    qc_path.append('extended')
if config['workflow']['reference_quant']:
    qc_path.append('reference')
if config['workflow']['exhaustive_ORFs']:
    qc_path.append('exhaustive')


rule all:
    input:
        ### Trimmed reads
        #expand("results/trimmed/{sample}_1.fq.gz", sample = SAMPLES),
        #expand("results/trimmed/{sample}_2.fq.gz", sample = SAMPLES),
        #multiqc_report = expand("results/{wf}/multiqc/multiqc_report.html", wf = qc_path)
        expand("results/{wf}.txt", wf = qc_path)
            

        


### snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
### snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg
### snakemake --rulegraph > rulegraph.dot && dot -Tsvg < rulegraph.dot > rulegraph.svg

### snakemake --use-conda --use-singularity -j 1 -r --verbose
### time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100


