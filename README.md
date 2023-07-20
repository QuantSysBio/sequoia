# Sequoia
#### Sequence Expression Quantification Unknown ORF discovery and isoform assembly


## Installation
```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install -c conda-forge singularity 
```

## Setup
1. Define common job parameters in `config.yml` and `features.yaml`
2. Place the paired reads in the `/data/reads_fq`. The first and second read files should end with `_1_fastq.gz` and `_2_fastq.gz` accordingly.  

## Execution
### Workflow check
`snakemake -j 1 -n`
### Workflow execution:
```
time snakemake --use-singularity --use-conda -j 1 --conda-frontend conda --resources load=100
```

## Slurm
Connect to the Mascot server or any calc node to submit the Slurm jobs. 

### Install conda in your home directory
Enter `bash Miniconda3-latest-Linux-x86_64.sh` and follow the instructions.
After that, create the SPIsnake environment as described under **Installation**.

### Clone repo + upload data
Enter `git clone https://github.com/QuantSysBio/Sequoia` to retrieve the latest code. You might need to [generate a token](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/) since GitHub recently removed the password authentication.
If necessary, deposit your data in the correct directory using `sftp`, `scp` or `rsync`. Instructions can be found in the [QSB getting started](https://pad.gwdg.de/s/JlkAOXJ2f#) document.

### Cluster execution
- Make sure you are in the correct directory  and on the correct node
- Sequoia is executed from a Bash screen session that prevents the job from crashing once you disconnect from `ssh`. Therefore, enter:
`screen -S sequoia`
- Activate the conda environment:
`conda activate snakemake`
- Submit the job to the `elbe` partition. (You can get an overview about which compute nodes are assigned to which partition by calling `sinfo`.)  
```
snakemake --use-singularity --cluster-config src/cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} --mem {cluster.mem} --job-name {cluster.job-name} -o {cluster.output} -D {cluster.chdir} --exclusive" --conda-frontend conda -j 3 -w 600 --restart-times 3 --resources load=100
```
- Detach from the screen session by pressing `Ctrl+a+d`. You can resume to the session to check the progress via `screen -r spisnake`

