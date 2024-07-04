
rule STAR_index:
    input:
        genome_fasta = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"]),
        reference_gtf = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"])
    output:
        directory(STAR_index_dir)
    benchmark: 
        "results/extended/benchmarks/STAR_index.txt"
    log: 
        "results/extended/logs/STAR_index.txt"
    conda: 
        "STAR.yaml"
    params: 
        n = config['program_specific']["star_cpus"]
    resources: 
        load = 100
    shell: 
        "mkdir -p {output} ; \
            STAR \
            --runThreadN {params.n} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --outFilterType BySJout \
            --alignSJoverhangMin 8 \
            --sjdbGTFfile {input.reference_gtf} \
            --genomeFastaFiles {input.genome_fasta} 2> {log}"


rule STAR_GTF_pass_1:
     input: 
        read1 = join(reads_trimmed, PATTERN_R1_trimmed),
        read2 = join(reads_trimmed, PATTERN_R2_trimmed),
        index = STAR_index_dir,
        reference_gtf = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"])
     output: 
        sj = 'results/extended/STAR_1st_pass/{sample}.SJ.out.tab',
        sams = temp('results/extended/STAR_1st_pass/{sample}.Aligned.out.sam')
     benchmark: 
        "results/extended/benchmarks/{sample}.align.json"
     log: 
        "results/extended/logs/{sample}.align_1pass.txt"
     conda: 
        "STAR.yaml"
     resources: 
        load = 100
     params: 
        n = config['program_specific']["star_cpus"]
     shell: "ulimit -n 10000 && \
            STAR \
         --runThreadN {params.n} \
         --genomeDir {input.index} \
         --readFilesIn {input.read1} {input.read2} \
         --readFilesCommand pigz -dc \
         --outFileNamePrefix results/extended/STAR_1st_pass/{wildcards.sample}. \
         --outSJfilterReads Unique \
         --sjdbGTFfile {input.reference_gtf} \2 > {log}"


rule merge_splice_junctions:
    input:
        sj = expand('results/extended/STAR_1st_pass/{sample}.SJ.out.tab', sample = SAMPLES)
    output:
        sjs = join('results/extended/STAR_1st_pass/SJ.out.pass1_merged.tab')
    log:
        'results/extended/logs/merge_splice_junctions.log'
    shell:
        # keep splice junctions with at least 3 uniquely mapped fragments per sample.
        "cat {input.sj} | awk '$7 >= 3' | cut -f1-4 | sort -u > {output.sjs}"


rule STAR_GTF_pass_2:
     input: 
        read1 = join(reads_trimmed, PATTERN_R1_trimmed),
        read2 = join(reads_trimmed, PATTERN_R2_trimmed),
        index = STAR_index_dir,
        sjs = 'results/extended/STAR_1st_pass/SJ.out.pass1_merged.tab',
        reference_gtf = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"])
     output: 
        #bam_tr= "results/STAR_2nd_pass/{sample}.Aligned.toTranscriptome.out.bam",
        bam = temp("results/extended/STAR_2nd_pass/{sample}.Aligned.sortedByCoord.out.bam")
     benchmark: 
        "results/extended/benchmarks/{sample}.align.json"
     log: 
        "results/extended/logs/{sample}.align_2pass.txt"
     conda: 
        "STAR.yaml"
     resources: 
        load = 100
     params: 
        n = config['program_specific']["star_cpus"]
     shell: "ulimit -n 10000 && \
            STAR \
            --quantMode GeneCounts \
            --runThreadN {params.n} \
            --genomeDir {input.index} \
            --readFilesIn {input.read1} {input.read2} \
            --outFileNamePrefix results/extended/STAR_2nd_pass/{wildcards.sample}. \
            --readFilesCommand pigz -dc \
            --outSAMattributes NH HI XS \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMattrIHstart 0 \
            --outSAMmapqUnique 60 \
            --outSJfilterReads Unique \
            --outFilterIntronMotifs None \
            --outReadsUnmapped None \
            --outFilterType BySJout \
            --alignSJoverhangMin 10 \
            --outFilterMultimapNmax 1  \
            --clip5pNbases 6 6 \
            --outFilterMismatchNmax 10 \
            --outFilterMismatchNoverLmax 0.1 \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outSAMattrRGline ID:GRPundef \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 \
            --chimNonchimScoreDropMin 10 \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --alignInsertionFlush Right \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --sjdbFileChrStartEnd {input.sjs} \
            --sjdbGTFfile {input.reference_gtf} \2 > {log}"

# for Salmon transcriptome             --quantMode TranscriptomeSAM GeneCounts \
