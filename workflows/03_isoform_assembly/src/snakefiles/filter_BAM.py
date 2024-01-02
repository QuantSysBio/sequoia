rule filter_BAM:
    input: 
        bam="results/STAR_2nd_pass/{sample}.Aligned.sortedByCoord.out.bam"
    output: 
        protected("results/filter_BAM/{sample}.Aligned.trimmed.out.bam")
    benchmark: 
        "results/benchmarks/{sample}.filter_BAM.txt"
    log: 
        "results/logs/{sample}.filter_BAM.txt"
    conda: 
        "filter_BAM.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "samtools view -@ {params.n} -b -h -F 4 -F 256 -F 512 -q 30 {input.bam} > {output} 2> {log}"


rule BuildBamIndex:
    input: 
        "results/filter_BAM/{sample}.Aligned.trimmed.out.bam"
    output: 
        "results/filter_BAM/{sample}.Aligned.trimmed.out.bai"
    benchmark: 
        "results/benchmarks/{sample}.BuildBamIndex.txt"
    log: 
        "results/logs/{sample}.BuildBamIndex.txt"
    conda: 
        "filter_BAM.yaml"
    params: 
        n=config["max_cores"]    
    shell: 
        "picard \
        BuildBamIndex \
        INPUT={input} 2> {log}"
