rule salmon_decoys:
    input:
        features["reference"]["genome_fasta"]
    output: 
        "results/trascriptome_quant/salmon/decoys.txt"
    benchmark: 
        "results/benchmarks/salmon_decoys.txt"
    log: 
        "results/logs/salmon_decoys.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = 1,
        mem = config["min_mem"],
        time = config["max_time"],
        load = 1
    shell: 
        """
        grep '^>' <{input} | cut -d ' ' -f 1 > {output}
        sed -i -e 's/>//g' {output}
        sed -i '/ERCC-/d' {output}
        """


rule salmon_gentrome:
    input:
        genome_fasta = features["reference"]["genome_fasta"],
        transcriptome_fasta = features["reference"]["transcriptome_fasta"]
    output: 
        gentrome="results/trascriptome_quant/salmon/gentrome.fasta"
    benchmark: 
        "results/benchmarks/salmon_gentrome.txt"
    log: 
        "results/logs/salmon_gentrome.txt"
    conda: 
        "R_remove_duplicates.yaml"
    resources:
        ncpus = 1,
        mem = config["min_mem"],
        time = config["max_time"],
        load = 10
    script: 
        "salmon_gentrome_ERCC.R"


rule salmon_index:
    input: 
        transcripts_fa = features["reference"]["transcriptome_fasta"],
        decoys = "results/trascriptome_quant/salmon/decoys.txt",
        gentrome = "results/trascriptome_quant/salmon/gentrome.fasta"
    output: 
        directory("results/trascriptome_quant/salmon/transcriptome_index"),
    benchmark: 
        "results/benchmarks/salmon_index.txt"
    log: 
        "results/logs/salmon_index.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"],
        load = 50
    params: 
        n = config["max_cpus"]
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
        index="results/trascriptome_quant/salmon/transcriptome_index",
        transcripts_fa = features["reference"]["transcriptome_fasta"],
        read1 = join(reads_trimmed, PATTERN_R1_trimmed),
        read2 = join(reads_trimmed, PATTERN_R2_trimmed)
    output: 
        tr_quant_sa="results/trascriptome_quant/salmon/{sample}/quant.sf"
    benchmark: 
        "results/benchmarks/salmon_quant_{sample}.txt"
    log: 
        "results/logs/salmon_quant_{sample}.txt"
    conda: 
        "salmon.yaml"
    resources:
        ncpus = config["salmon_cpus"],
        mem = config["salmon_mem"],
        time = config["max_time"],
        load = 25
    params: 
        n = config["salmon_cpus"],
        salmon_Gibbs_samples = config["salmon_Gibbs_samples"]
    shell: 
        "salmon quant \
            --libType A \
            -i {input.index} \
            -p {params.n} \
            -1 {input.read1} \
            -2 {input.read2} \
            -o results/trascriptome_quant/salmon/{wildcards.sample}/ \
            --validateMappings \
            --mimicBT2 \
            --rangeFactorizationBins 4 \
            --seqBias \
            --gcBias \
            --reduceGCMemory \
            --posBias \
            --numGibbsSamples {params.salmon_Gibbs_samples} \
            2> {log}"
