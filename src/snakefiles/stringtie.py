# Define SAMPLES list somewhere in your config or Snakefile


# Rule to process GTF with StringTie
rule StringTie_GTF:
    input:
        bam = "results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam",
        reference_gtf = "data/references/gencode.v45.primary_assembly.annotation_copy.gtf" if config["reference"]['ERCC_normalization'] else config["reference"]["reference_gtf"]
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
        """
        stringtie \
            -G {input.reference_gtf} \
            {input.bam} \
            -p {params.n} \
            -o {output} \
            --conservative \
            -m 200 \
            -v 2> {log}
        """

# Rule to merge GTF files
rule merge_GTF:
    input:
        expand("results/extended/stringtie_assemble/{sample}-stringtie.gtf", sample=SAMPLES),
        reference_gtf = "data/references/gencode.v45.primary_assembly.annotation_copy.gtf" if config["reference"]['ERCC_normalization'] else config["reference"]["reference_gtf"]
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
        """
        stringtie \
            --merge \
            -G {input.reference_gtf} \
            -o {output.merged_gtf} \
            -p {params.n} \
            -c 0 \
            -m 50 \
            -T 0 \
            -f .05 \
            {input} 2> {log}
        """

# Rule for gffcompare
rule gff_compare:
    input:
        reference_gtf = "data/references/gencode.v45.primary_assembly.annotation_copy.gtf" if config["reference"]['ERCC_normalization'] else config["reference"]["reference_gtf"],
        genome_fasta = "data/references/GRCh38.primary_assembly.genome_copy.fa" if config["reference"]['ERCC_normalization'] else config["reference"]["genome_fasta"],
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
        """
        gffcompare \
            -r {input.reference_gtf} \
            -s {input.genome_fasta} \
            -D \
            -V \
            -A \
            -X \
            -K \
            -o results/extended/stringtie_assemble/gffcompare/gffcmp \
            {input.merged_gtf} 2> {log}
        """

# Rule for StringTie quantification
rule stringtie_quant:
    input:
        bam = "results/extended/filter_BAM/{sample}.Aligned.trimmed.out.bam",
        gff_compare = "results/extended/stringtie_assemble/gffcompare/gffcmp.combined.gtf"
    output:
        tr_quant_st = "results/extended/trascriptome_quant/stringtie/{sample}/{sample}_t_data.ctab"
    benchmark:
        "results/extended/benchmarks/stringtie_quant_{sample}.txt"
    log:
        "results/extended/logs/stringtie_quant_{sample}.txt"
    conda:
        "stringtie.yaml"
    params:
        n = config['program_specific']["stringtie_cpus"]
    shell:
        """
        stringtie \
            -G {input.gff_compare} \
            {input.bam} \
            -p {params.n} \
            -o {output.tr_quant_st} \
            -e -B \
            -v 2> {log}
        """
