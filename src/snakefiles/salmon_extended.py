rule salmon_decoys_exttended:
    input:
        branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"])
    output: 
        "results/extended/trascriptome_quant/salmon/decoys.txt"
    benchmark: 
        "results/extended/benchmarks/salmon_decoys.txt"
    log: 
        "results/extended/logs/salmon_decoys.txt"
    conda: 
        "salmon.yaml"
    shell: 
        """
        grep '^>' <{input} | cut -d ' ' -f 1 > {output}
        sed -i -e 's/>//g' {output}
        sed -i '/ERCC-/d' {output}
        """


rule salmon_gentrome_extended:
    input:
        genome_fasta=branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"]),
        transcriptome_fasta = "results/extended/tr_2_prot/merged_reference.transcripts.fasta"
    output: 
        gentrome="results/extended/trascriptome_quant/salmon/gentrome.fasta"
    benchmark: 
        "results/extended/benchmarks/salmon_gentrome.txt"
    log: 
        "results/extended/logs/salmon_gentrome.txt"
    conda: 
        "R_remove_duplicates.yaml"
    resources: 
        load = 50
    script: 
        "salmon_gentrome_ERCC.R"


rule salmon_index_extended:
    input: 
        transcripts_fa = "results/extended/tr_2_prot/merged_reference.transcripts.fasta",
        decoys = "results/extended/trascriptome_quant/salmon/decoys.txt",
        gentrome = "results/extended/trascriptome_quant/salmon/gentrome.fasta"
    output: 
        directory("results/extended/trascriptome_quant/salmon/transcriptome_index"),
    benchmark: 
        "results/extended/benchmarks/salmon_index.txt"
    log: 
        "results/extended/logs/salmon_index.txt"
    conda: 
        "salmon.yaml"
    resources: 
        load = 100
    params: 
        n=config['program_specific']["salmon_cpus"]
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
rule salmon_quant_extended:
    input: 
        index="results/extended/trascriptome_quant/salmon/transcriptome_index",
        transcripts_fa="results/extended/tr_2_prot/merged_reference.transcripts.fasta",
        read1 = join(reads_trimmed, PATTERN_R1_trimmed),
        read2 = join(reads_trimmed, PATTERN_R2_trimmed)
    output: 
        tr_quant_sa="results/extended/trascriptome_quant/salmon/{sample}/quant.sf"
    benchmark: 
        "results/extended/benchmarks/salmon_quant_{sample}.txt"
    log: 
        "results/extended/logs/salmon_quant_{sample}.txt"
    conda: 
        "salmon.yaml"
    resources: 
        load = 100
    params: 
        n=config['program_specific']["salmon_cpus"]
    shell: 
        "salmon quant \
            --libType A \
            -i {input.index} \
            -p {params.n} \
            -1 {input.read1} \
            -2 {input.read2} \
            -o results/extended/trascriptome_quant/salmon/{wildcards.sample}/ \
            --validateMappings \
            --mimicBT2 \
            --rangeFactorizationBins 4 \
            --seqBias \
            --gcBias \
            --reduceGCMemory \
            --posBias \
            --numGibbsSamples 1000 \
            2> {log}"


### Alternative: Alignment-based quantifiction
rule salmon_quant_BAM:
    input: 
        transcripts_fa="results/extended/tr_2_prot/merged_reference.transcripts.fasta",
        reference_gtf="results/extended/tr_2_prot/merged_reference.transcripts.gff",
        bam="results/extended/STAR_2nd_pass/{sample}.Aligned.toTranscriptome.out.bam"
    output: 
        tr_quant_BAM="results/extended/trascriptome_quant/salmon_BAM/{sample}/quant.sf"
    benchmark: 
        "results/extended/benchmarks/salmon_quant_BAM_{sample}.txt"
    log: 
        "results/extended/logs/salmon_quant_BAM_{sample}.txt"
    conda: 
        "salmon.yaml"
    resources: 
        load = 100
    params: 
        n=config['program_specific']["salmon_cpus"]
    shell: 
        "salmon quant \
            --libType A \
            -p {params.n} \
            -t {input.transcripts_fa} \
            -g {input.reference_gtf} \
            -a {input.bam} \
            -o results/extended/trascriptome_quant/salmon_BAM/{wildcards.sample}/ \
            --seqBias \
            --gcBias \
            --posBias \
            --reduceGCMemory \
            --numGibbsSamples 1000 \
            --gencode \
            2> {log}"
