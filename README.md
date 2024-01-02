# Sequoia
#### Sequence Expression Quantification Unknown ORF discovery and Isoform Assembly
## Description
Sequoia allows for mass spectrometry search space definition and facilitates multi-omics integration. Three workflows are available: 
- (1) Reference transcriptome quantification from RNA-seq data
- (2) Reference-guided transcriptome assembly followed by ORF search on the new transcripts and quantification of the augmented transcriptome. 
- (3) Exhaustive ORF definition can be done with genome and transcriptome information.

## Installation
```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge singularity 
```

## Setup
1. Define common job parameters in `config.yml` and `features.yaml`
2. Place the paired reads in the `/data/reads_fq`. The first and second read files should end with `_1_fastq.gz` and `_2_fastq.gz` accordingly.
3. Download the [GENCODE annotation](https://www.gencodegenes.org/) to `data/reference` and match the files to `features.yaml`
4. (for isoform assembly workflow) download the [Pfam data](https://www.ebi.ac.uk/interpro/download/pfam/) to `data/hmm`
5. (for exhaustive ORF workflow) the `00_BSgenome_Gencode.R` script allows to use the external genomes starting from a fasta file and demonstrates the setup for GENCODE. To skip this step, use the [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html) or any other genome as a Biostrings object. 
## Execution
### Workflow check
`snakemake -j 1 -n`
### Workflow execution:
```
time snakemake --use-singularity --use-conda -j 1 --conda-frontend conda --resources load=100
```

### Install conda in your home directory
Enter `bash Miniconda3-latest-Linux-x86_64.sh` and follow the instructions.
After that, create the SPIsnake environment as described under **Installation**.

### Clone repo + upload data
Enter `git clone https://github.com/QuantSysBio/sequoia` to retrieve the latest code. You might need to [generate a token](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/) since GitHub recently removed the password authentication.

### Cluster execution with Slurm
- Sequoia is executed from a Bash screen session that prevents the job from terminating once you disconnect from `ssh`. Therefore, enter:
`screen -S sequoia`
- Activate the conda environment:
`conda activate snakemake`
- fill in the cluster configuration:`src/cluster.yaml`
```
snakemake --use-singularity --use-conda --cluster-config src/cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} --mem {cluster.mem} --job-name {cluster.job-name} -o {cluster.output} -D {cluster.chdir} --exclusive" --conda-frontend conda -j 1 -w 60 --resources load=100
```
- Detach from the screen session by pressing `Ctrl+a+d`. You can resume to the session to check the progress via `screen -r sequoia`

