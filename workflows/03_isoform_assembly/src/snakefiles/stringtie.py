rule StringTie_GTF:
        input: 
            bam="results/filter_BAM/{sample}.Aligned.trimmed.out.bam",
            reference_gtf=features["reference"]["reference_gtf"]
        output: 
            "results/stringtie_assemble/{sample}-stringtie.gtf"
        benchmark: 
            "results/benchmarks/{sample}.StringTie.txt"
        log: 
            "results/logs/{sample}.StringTie.txt"
        conda: 
            "stringtie.yaml"
        params: 
            n=config["max_cores"]
        resources: 
            load = 20
        shell: "stringtie \
                -G {input.reference_gtf} \
                {input.bam} \
                -p {params.n} \
                -o {output} \
                --conservative \
                -m 200 \
                -v 2> {log}"


rule merge_GTF:
    input: 
        expand("results/stringtie_assemble/{sample}-stringtie.gtf", sample=SAMPLES),
        reference_gtf=features["reference"]["reference_gtf"]
    output: 
        merged_gtf="results/stringtie_assemble/merged.gtf"
    benchmark: 
        "results/benchmarks/merge.txt"
    log: 
        "results/logs/merge.txt"
    conda: 
        "stringtie.yaml"
    params: 
        n=config["max_cores"]
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
    input: reference_gtf=features["reference"]["reference_gtf"],
           genome_fasta=features["reference"]["genome_fasta"],
           merged_gtf="results/stringtie_assemble/merged.gtf"
    output: 
        gff_compare="results/stringtie_assemble/gffcompare/gffcmp.combined.gtf",
        gffcomp_tracking = "results/stringtie_assemble/gffcompare/gffcmp.tracking"
    benchmark: 
        "results/benchmarks/gff_compare.txt"
    log: 
        "results/logs/gff_compare.txt"
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
            -o results/stringtie_assemble/gffcompare/gffcmp \
            {input.merged_gtf}"


rule stringtie_quant:
    input: 
        bam="results/filter_BAM/{sample}.Aligned.trimmed.out.bam",
        gff_compare="results/stringtie_assemble/gffcompare/gffcmp.combined.gtf"
    output: 
        tr_quant_st="results/trascriptome_quant/stringtie/{sample}/{sample}_t_data.ctab"
    benchmark: 
        "results/benchmarks/stringtie_quant_{sample}.txt"
    log: 
        "results/logs/stringtie_quant_{sample}.txt"
    conda: 
        "stringtie.yaml"
    params: 
        n=config["max_cores"]
    shell: 
        "stringtie  \
        -G {input.gff_compare} \
        {input.bam} \
        -p {params.n} \
        -o {output} \
        -e -B \
        -v 2> {log}"