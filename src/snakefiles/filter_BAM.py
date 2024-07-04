# Rule to filter BAM files using samtools
rule filter_BAM:
    input:
        bam = "results/extended/STAR_2nd_pass/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        protected("results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam")
    benchmark:
        "results/extended/benchmarks/{sample}.filter_BAM.txt"
    log:
        "results/extended/logs/{sample}.filter_BAM.txt"
    conda:
        "filter_BAM.yaml"
    params:
        n = config['resources']["max_cpus"]
    shell:
        """
        samtools view -@ {params.n} -b -h -F 4 -F 256 -F 512 -q 30 {input.bam} > {output} 2> {log}
        """

# Rule to build BAM index using Picard
rule BuildBamIndex:
    input:
        bam = "results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam"
    output:
        bai = "results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bai"
    benchmark:
        "results/extended/benchmarks/{sample}.BuildBamIndex.txt"
    log:
        "results/extended/logs/{sample}.BuildBamIndex.txt"
    conda:
        "filter_BAM.yaml"
    params:
        n = config['resources']["max_cpus"]
    shell:
        """
        picard BuildBamIndex \
        INPUT={input.bam} \
        OUTPUT={output.bai} \
        2> {log}
        """
