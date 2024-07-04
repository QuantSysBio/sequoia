rule StringTie_GTF:
        input: 
            bam="results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam",
            reference_gtf = "data/references/gencode.v45.primary_assembly.annotation_copy.gtf" if config["reference"]['include_ERCC'] else config["reference"]["reference_gtf"]
        output:
            "results/extended/stringtie_assemble/{sample}-stringtie.gtf"
        benchmark:
            "results/extended/benchmarks/{sample}.StringTie.txt"
        log:
            "results/extended/logs/{sample}.StringTie.txt"
        conda: 
            "stringtie.yaml"
        params: 
            n = config['program_specific']["stringtie_cpus"]
        resources: 
            load = 20
        shell:
            "stringtie \
            -G {input.reference_gtf} \
            {input.bam} \
            -p {params.n} \
            -o {output} \
            --conservative \
            -m 200 \
            -v 2> {log}
            "

rule merge_GTF:
    input:
        expand("results/extended/stringtie_assemble/{sample}-stringtie.gtf", sample=SAMPLES),
        reference_gtf = "data/references/gencode.v45.primary_assembly.annotation_copy.gtf" if config["reference"]['include_ERCC'] else config["reference"]["reference_gtf"]
    output:
        merged_gtf = "results/extended/stringtie_assemble/merged.gtf"
    benchmark:
        "results/extended/benchmarks/merge.txt"
    log:
        "results/extended/logs/merge.txt"
    conda:
        "stringtie.yaml"
    params:
        n = config['program_specific']["stringtie_cpus"]
    resources:
            load = 20
    shell:
        "stringtie \
            --merge \
            -G {input.reference_gtf} \
            -o {output} \
            -p {params.n} \
            -c 0 \
            -m 50 \
            -T 0 \
            -f .05 \
            {input} 2> {log}"


rule gff_compare:
    input: reference_gtf = branch(config["reference"]['include_ERCC'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"]),
           genome_fasta = branch(config["reference"]['include_ERCC'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"]),
           merged_gtf = "results/extended/stringtie_assemble/merged.gtf"
    output: 
        gff_compare = "results/extended/stringtie_assemble/gffcompare/gffcmp.combined.gtf",
        gffcomp_tracking = "results/extended/stringtie_assemble/gffcompare/gffcmp.tracking"
    benchmark: 
        "results/extended/benchmarks/gff_compare.txt"
    log: 
        "results/extended/logs/gff_compare.txt"
    conda: 
        "stringtie.yaml"
    resources: 
        load = 20
    shell: 
        "gffcompare \
            -r {input.reference_gtf} \
            -s {input.genome_fasta} \
            -D \
            -V \
            -A \
            -X \
            -K \
            -o results/extended/stringtie_assemble/gffcompare/gffcmp \
            {input.merged_gtf}"


rule stringtie_quant:
    input: 
        bam="results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam",
        gff_compare="results/extended/stringtie_assemble/gffcompare/gffcmp.combined.gtf"
    output: 
        tr_quant_st="results/extended/trascriptome_quant/stringtie/{sample}/{sample}_t_data.ctab"
    benchmark: 
        "results/extended/benchmarks/stringtie_quant_{sample}.txt"
    log: 
        "results/extended/logs/stringtie_quant_{sample}.txt"
    conda: 
        "stringtie.yaml"
    params: 
        n = config['program_specific']["stringtie_cpus"]
    shell: 
        "stringtie  \
        -G {input.gff_compare} \
        {input.bam} \
        -p {params.n} \
        -o {output} \
        -e -B \
        -v 2> {log}"