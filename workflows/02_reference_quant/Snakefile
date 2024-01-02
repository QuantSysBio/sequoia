shell.prefix("set -euo pipefail;")
shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

container: "docker://snakemake/snakemake:v6.8.0"
min_version("6.0")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
include: snakefiles + "folders.py"
include: snakefiles + "reads_wildcards.py"
include: snakefiles + "pre_process.py"
include: snakefiles + "salmon.py"
include: snakefiles + "transcriptome_2_proteome.py"
include: snakefiles + "multiqc.py"


rule all:
    input:
        ### Trimmed reads
        expand("results/trimmed/{sample}_1.fq.gz", sample = SAMPLES),
        expand("results/trimmed/{sample}_2.fq.gz", sample = SAMPLES),

        ### Expressed proteome
        expr_prot="results/tr_2_prot/proteome_expressed.fasta",
        
        ### Salmon: selective alignment
        tr_quant_sa=expand("results/trascriptome_quant/salmon/{sample}/quant.sf", sample = SAMPLES),

        ### MultiQC
        multiqc_report="results/multiqc/multiqc_report.html"
        


### snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
### snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg
### snakemake --rulegraph > rulegraph.dot && dot -Tsvg < rulegraph.dot > rulegraph.svg

### snakemake --use-conda --use-singularity -j 1 -r --verbose
### time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100


