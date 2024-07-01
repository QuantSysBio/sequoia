# Sequoia
#### Sequence Expression Quantification Unknown ORF discovery and Isoform Assembly

Sequoia is a pipeline for mass spectrometry search space definition and facilitates multi-omics integration. Three workflows are available:  
(1) Reference transcriptome quantification from RNA-seq data.  
(2) Reference-guided transcriptome assembly followed by ORF prediction on the new transcripts and quantification of the augmented transcriptome.  
(3) Exhaustive ORF search with genome and transcriptome annotation.  

## Installation
```
conda create -c conda-forge -c bioconda -n snakemake snakemake==8.10.4
conda activate snakemake
conda install -c conda-forge singularity==3.8.6
```

## Setup
1. Define job parameters in `config.yml`
2. Place the paired reads in the `/data/reads_fq`. The first and second read files should end with `_1.fq.gz` and `_2.fq.gz` accordingly.  
3. (optional) Fill in the Slurm cluster configuration:`src/cluster.yaml`

## Execution
### Workflow check
```
snakemake -j 1 -n
```

### Workflow execution:
```
time snakemake --use-singularity --use-conda -j 1 --conda-frontend mamba --resources load=100
```

### Cluster execution with Slurm
- Sequoia can be executed from a Bash screen session that prevents the job from terminating once you disconnect from `ssh`. Therefore, enter:
  `screen -S sequoia`
- Activate the conda environment:
  `conda activate snakemake`
- fill in the cluster configuration:`src/cluster.yaml`
```
snakemake --use-singularity --use-conda --cluster-config src/cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} --mem {cluster.mem} --job-name {cluster.job-name} -o {cluster.output} -D {cluster.chdir} --exclusive" --conda-frontend conda -j 1 -w 60 --resources load=100
```
- Detach from the screen session by pressing `Ctrl+a+d`. You can resume to the session to check the progress via `screen -r sequoia`

### Input description
| Parameter    | Meaning    | 
  | ----------- | ----------- | 
  | workflow | Whether to execute a given workflow: True or False | 
  | reference | Path to the reference files | 
  | resources | General resources to run tools if their specific resources are not set in program_specific section | 
  | slurm | Specific settings for running the pipeline on slurm |
  | program_specific | Set parameters specific for each tool separately |

### All options
  
  | Parameter    | Meaning    | type |
  | -----------  | -----------| -----|
  | max_mem | Maximum memory for all tools unless specified otherwise | Int |
| min_mem | Minimum memory for all tools unless specified otherwise | Int |
  | min_cpus | Minimum CPUs for all tools unless specified otherwise | Int |
| max_cpus | Maximum CPUs for all tools unless specified otherwise | Int |
  | exhaustive_ORFs | Enable the workflow to search for all possible ORFs across transcriptome and genome | Boolean |
  | reference_quant | Enable the workflow to quantify the reference transcriptome using <a href=https://doi.org/10.1038/nmeth.4197>`Salmon`</a> | Boolean |
  | isoform_assembly | Enable the workflow to align RNA-seqs to reference by <a href=https://doi.org/10.1093/bioinformatics/bts635>`STAR`</a> and assemble Reference-guided transcriptome using <a href=https://ccb.jhu.edu/software/stringtie>`Stringtie`</a> followed by ORF prediction by <a href=https://github.com/TransDecoder/TransDecoder>`TransDecoder`</a> on the new transcripts and quantify the augmented transcriptome by `Salmon`| Boolean |
  | genome_fasta | path to the reference genome | String |
  | transcriptome_fasta | path to the reference transcriptome | String |
  | gencode_translation_sequences_fasta | path to the reference proteome | String |
  | reference_gtf | path to the reference gtf annotation | String |
  | reference_gff3 | path to the reference gff3 annotation | String |
  | gencode_metadata_TrEMBL |path to the reference mapping to TrEMBL | String |
  | gencode_metadata_SwissProt |path to the reference  mapping to SwissProt | String |
  | ERCC_normalization | Whether to normalize the quantification using <a href=https://www.nist.gov/programs-projects/external-rna-controls-consortium> external ERCC</a> spike in sequences | String |
  | external_sequences | path to the external sequences to be included into the quantification  | String |
  | fastp_cpus | Maximum CPUs for <a href=https://doi.org/10.1002/imt2.107>`Fastp tool`</a> | Int |
  | star_cpus | Maximum CPUs for <a herf=https://doi.org/10.1093/bioinformatics/bts635>STAR tool</a> | Int |
  | stringtie_cpus | Maximum CPUs for Stringtie tool | Int |
  | salmon_cpus | Maximum CPUs for Salmon tool | Int |
  | salmon_mem | Maximum memory for Salmon tool | Numeric |
  | salmon_Gibbs_samples | Number of samples to be drawn from the posterior by Salmon tool | Int |
  | min_filt_samples | Minimal number of bio.reps required for transcript identification | Int |
  | min_filt_counts | Minimal number of read counts required for transcript identification  | Numeric |
  | phred_score | Minimum quality score for trimming raw reads | Numeric |
  | min_orf_length | Minimum length of an ORF to be considered in exhaustive workflow  | Int |
  | start_codons | Which start codons to consider during the exhaustive ORF search (separated by only a comma like: ATG,CTG,GTG,ATC,ACG as set by default in config file)  | String |
  

### Output description
Results of each workflow is saved in results folder under a folder named according to the workflow name (/results/('Exhaustive', 'refrence' and 'extended')). Moreover, final quantification results and QC output can be seen under the 'transcriptome_quant' and 'multiqc' folders, respectively, which all located under the mentioned results folder for each workflow. Moreover, benchmarks and logs folder contained information regarding tools performance and logs created by them, respectively. 


```bash

├── extended
│   ├── benchmarks
│   ├── filter_BAM
│   ├── logs
│   ├── multiqc
│   │   ├── multiqc_data
│   │   └── multiqc_report.html
│   ├── STAR_1st_pass
│   ├── STAR_2nd_pass
│   ├── STAR_index
│   ├── stringtie_assemble
│   ├── tr_2_prot
│   └── trascriptome_quant
│       ├── gene_expression.csv
│       ├── gene_expression_tpm.csv
│       ├── transcript_expression.csv
│       └── transcript_expression_tpm.csv
├── extended.txt (Indicator of isoform_assembly workflow is finished)
├── logs_fastp
├── reference.txt (Indicator of reference_quant workflow is finished)
├── refrence
│   ├── benchmarks
│   ├── logs
│   ├── multiqc
│   │   └── multiqc_report.html
│   ├── tr_2_prot
│   └── transcriptome_quant
│       ├── gene_expression.csv
│       ├── gene_expression_tpm.csv
│       ├── salmon
│       ├── transcript_expression.csv
│       └── transcript_expression_tpm.csv
├── Exhaustive
│   ├── benchmarks
│   ├── CDS
│   ├── cryptic
│   └── logs
├── exhaustive.txt (Indicator of exhaustive_ORFs workflow is finished)
└── trimmed (Output trimmed sequences of fastp tool saved here)

 ```
