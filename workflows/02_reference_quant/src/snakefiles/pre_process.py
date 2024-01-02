rule fastp_pe:
    input:
        read1 = join(reads_fq, PATTERN_R1),
        read2 = join(reads_fq, PATTERN_R2)
    output:
        read1 = "results/trimmed/{sample}_1.fq.gz",
        read2 = "results/trimmed/{sample}_2.fq.gz"
    params: 
        n = config["fastp_cpus"],
        phred = config["phred_score"]
    conda: 
        "pre_process.yaml"
    resources:
        ncpus = config["fastp_cpus"],
        mem = config["min_mem"],
        time = config["max_time"],
        load = 5
    log: 
        json = "results/logs/{sample}.fastp.json"
    shell:
        "fastp -i {input.read1} \
        -I {input.read2} \
        -o {output.read1} \
        -O {output.read2} \
        -j {log.json} \
        --qualified_quality_phred {params.phred} \
        --correction \
        --verbose \
        --low_complexity_filter \
        --thread {params.n}"
