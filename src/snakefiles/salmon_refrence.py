          
rule salmon_decoys:
    input:
        branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"])
    output: 
        "results/refrence/transcriptome_quant/salmon/decoys.txt"
    benchmark: 
        "results/refrence/benchmarks/salmon_decoys.txt"
    log: 
        "results/refrence/logs/salmon_decoys.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = 1,
        mem = config["resources"]["min_mem"],
        load = 1
    shell: 
        """
        grep '^>' <{input} | cut -d ' ' -f 1 > {output}
        sed -i -e 's/>//g' {output}
        sed -i '/ERCC-/d' {output}
        """


rule salmon_gentrome:
    input:
        genome_fasta = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"]),
        transcriptome_fasta = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.transcripts_copy.fa',
                              otherwise = config["reference"]["transcriptome_fasta"])
    output: 
        gentrome="results/refrence/transcriptome_quant/salmon/gentrome.fasta"
    benchmark: 
        "results/refrence/benchmarks/salmon_gentrome.txt"
    log: 
        "results/refrence/logs/salmon_gentrome.txt"
    conda: 
        "R_remove_duplicates.yaml"
    resources:
        ncpus = 1,
        mem = config["resources"]["min_mem"],
        load = 10
    script: 
        "salmon_gentrome_ERCC.R"


rule salmon_index:
    input: 
        transcripts_fa = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.transcripts_copy.fa',
                              otherwise = config["reference"]["transcriptome_fasta"]),
        decoys = "results/refrence/transcriptome_quant/salmon/decoys.txt",
        gentrome = "results/refrence/transcriptome_quant/salmon/gentrome.fasta"
    output: 
        directory("results/refrence/transcriptome_quant/salmon/transcriptome_index"),
    benchmark: 
        "results/refrence/benchmarks/salmon_index.txt"
    log: 
        "results/refrence/logs/salmon_index.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = config["resources"]["max_cpus"],
        mem = config["resources"]["max_mem"],
        load = 50
    params: 
        n = config["resources"]["max_cpus"]
    shell: 
        "salmon index \
        --transcripts {input.gentrome} \
        --index {output} \
        --decoys {input.decoys} \
        --kmerLen 31 \
        -p {params.n}\
        --keepDuplicates \
        --gencode \
        2> {log}"


### Selective alignment
rule salmon_quant:
    input: 
        index = "results/refrence/transcriptome_quant/salmon/transcriptome_index",
        transcripts_fa = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gff3',
                              otherwise = config["reference"]["transcriptome_fasta"]),
        read1 = join(reads_trimmed, PATTERN_R1_trimmed),
        read2 = join(reads_trimmed, PATTERN_R2_trimmed)
    output: 
        tr_quant_sa = "results/refrence/transcriptome_quant/salmon/{sample}/quant.sf"
    benchmark: 
        "results/refrence/benchmarks/salmon_quant_{sample}.txt"
    log: 
        "results/refrence/logs/salmon_quant_{sample}.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = config["program_specific"]["salmon_cpus"],
        mem = config["program_specific"]["salmon_mem"],
        # load = 25
    params: 
        n = config["program_specific"]["salmon_cpus"],
        salmon_Gibbs_samples = config["program_specific"]["salmon_Gibbs_samples"]
    shell: 
        "salmon quant \
        --libType A \
        -i {input.index} \
        -p {params.n} \
        -1 {input.read1} \
        -2 {input.read2} \
        -o results/refrence/transcriptome_quant/salmon/{wildcards.sample}/ \
        --validateMappings \
        --mimicBT2 \
        --rangeFactorizationBins 4 \
        --seqBias \
        --gcBias \
        --reduceGCMemory \
        --posBias \
        --numGibbsSamples {params.salmon_Gibbs_samples} \
        2> {log}"
