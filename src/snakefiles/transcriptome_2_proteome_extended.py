rule transcript_sequences:
    """
    Cluster gffcompare outputs and get transcript sequences
    """
    input: 
        gff_compare="results/extended/stringtie_assemble/gffcompare/gffcmp.combined.gtf",
        genome_fasta = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"]),
        merged_gtf="results/extended/stringtie_assemble/merged.gtf"
    output: 
        mer_fa="results/extended/tr_2_prot/gffcmp.transcripts.fasta",
        mer_gtf="results/extended/tr_2_prot/clust_reference.transcripts.gtf"
    benchmark: 
        "results/extended/benchmarks/transcript_sequences.json"
    log: 
        "results/extended/logs/transcript_sequences.txt"
    conda: 
        "transcriptome_2_proteome.yaml"
    params: 
        n=config['resources']["max_cpus"]
    shell: 
        "gffread {input.gff_compare} \
        -g {input.genome_fasta} \
        -F \
        -P --adj-stop \
        --no-pseudo \
        --force-exons \
        --keep-genes \
        --keep-comments \
        --merge -K -Q -Y \
        --t-adopt \
        -v -E \
        -T -o {output.mer_gtf} \
        -w {output.mer_fa}  2> {log}"


rule remove_duplicate_transcripts:
    """
    Create a hybrid transcriptome annotation from reference and stringtie output
    """
    input:
        mer_gtf="results/extended/tr_2_prot/clust_reference.transcripts.gtf",
        mer_fa="results/extended/tr_2_prot/gffcmp.transcripts.fasta",
        ref_gtf=branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"]),        
        ref_fa=branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.transcripts_copy.fa',
                              otherwise = config["reference"]["transcriptome_fasta"])
    output:
        fasta="results/extended/tr_2_prot/merged_reference.transcripts.fasta",
        denovo_fasta = "results/extended/tr_2_prot/denovo.transcripts.fasta",
        gff_final="results/extended/tr_2_prot/merged_reference.transcripts.gff",
        gtf_denovo = "results/extended/tr_2_prot/denovo.transcripts.gtf",
    benchmark: 
        "results/extended/benchmarks/remove_duplicate_transcripts.json"
    log: 
        "results/extended/logs/remove_duplicate_transcripts.txt"
    conda: 
        "R_remove_duplicates.yaml"
    script: 
        "remove_duplicate_transcripts.R"


rule transdecoder_LongOrfs:
    """
    Find longest ORFs on novel transcripts
    """
    input:
        fasta = "results/extended/tr_2_prot/denovo.transcripts.fasta"    
    output:
        "results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    benchmark: 
        "results/extended/benchmarks/transdecoder_LongOrfs.json"
    log: 
        "results/extended/logs/transdecoder_LongOrfs.txt"
    conda: 
        "transdecoder.yaml"
    shell:
        "TransDecoder.LongOrfs \
        -t {input.fasta} \
        -O results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/ \
        -m 100 &> {log}"


rule transdecoder_copyfiles:
    input: 
        pep = "results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: 
        pep="results/extended/tr_2_prot/longest_orfs.pep",
        gff3="results/extended/tr_2_prot/longest_orfs.gff3",
        cds="results/extended/tr_2_prot/longest_orfs.cds",
        dat="results/extended/tr_2_prot/base_freqs.dat"
    benchmark: 
        "results/extended/benchmarks/transdecoder_copyfiles.json"
    log: 
        "results/extended/logs/transdecoder_copyfiles.txt"
    params:
        gff3 = "results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/longest_orfs.gff3",
        cds = "results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/longest_orfs.cds",
        dat = "results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/base_freqs.dat"
    shell:
        "cp -f {input.pep} {output.pep} && \
        cp -f {params.gff3} {output.gff3} && \
        cp -f {params.cds} {output.cds} && \
        cp -f {params.dat} {output.dat}"

rule download_hmm:
    input:
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    output:
        "data/hmm/Pfam-A.hmm.gz"
    shell:
        "wget -O {output} https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

rule decompress_hmm:
    input:
        "data/hmm/Pfam-A.hmm.gz"
    output:
        "data/hmm/Pfam-A.hmm"
    shell:
        "gzip -dk {input}"

rule hmmpress:
    input: 
        "data/hmm/Pfam-A.hmm"
    output: 
        h3f="data/hmm/Pfam-A.hmm.h3f",
        h3i="data/hmm/Pfam-A.hmm.h3i",
        h3m="data/hmm/Pfam-A.hmm.h3m",
        h3p="data/hmm/Pfam-A.hmm.h3p"
    benchmark: 
        "results/extended/benchmarks/hmmpress.json"
    log: 
        "results/extended/logs/hmmpress.txt"
    conda: 
        "hmmscan.yaml"
    resources: 
        load = 50
    params: 
        h3f="results/tr_2_prot/Pfam-A.hmm.h3f",
        h3i="results/tr_2_prot/Pfam-A.hmm.h3i",
        h3m="results/tr_2_prot/Pfam-A.hmm.h3m",
        h3p="results/tr_2_prot/Pfam-A.hmm.h3p",
        hmm="results/tr_2_prot/Pfam-A.hmm"
    shell: 
        "hmmpress {input} \
        2> {log}"


rule hmm_db_copyfiles:
    input: 
        h3p="data/hmm/Pfam-A.hmm.h3p"
    output: 
        h3f="results/extended/tr_2_prot/Pfam-A.hmm.h3f",
        h3i="results/extended/tr_2_prot/Pfam-A.hmm.h3i",
        h3m="results/extended/tr_2_prot/Pfam-A.hmm.h3m",
        h3p="results/extended/tr_2_prot/Pfam-A.hmm.h3p",
        hmm="results/extended/tr_2_prot/Pfam-A.hmm"
    benchmark: 
        "results/extended/benchmarks/hmm_db_copyfiles.json"
    log: 
        "results/extended/logs/hmm_db_copyfiles.txt"
    params:
        h3f="data/hmm/Pfam-A.hmm.h3f",
        h3i="data/hmm/Pfam-A.hmm.h3i",
        h3m="data/hmm/Pfam-A.hmm.h3m",
        h3p="data/hmm/Pfam-A.hmm.h3p",
        hmm="data/hmm/Pfam-A.hmm"
    shell:
        "cp -f {input.h3p} {output.h3p} && \
        cp -f {params.h3f} {output.h3f} && \
        cp -f {params.h3i} {output.h3i} && \
        cp -f {params.h3m} {output.h3m} && \
        cp -f {params.hmm} {output.hmm}"


# Can be executed with hmmscan. See:
# http://cryptogenomicon.org/hmmscan-vs-hmmsearch-speed-the-numerology.html
rule hmmscan:
    """
    Look for domains with Hmmer
    """
    input: 
        longest_orfs="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        hmm="results/extended/tr_2_prot/Pfam-A.hmm"
    output: 
        "results/extended/tr_2_prot/hmmscan/pfam.domtblout"
    benchmark: 
        "results/extended/benchmarks/hmmscan.json"
    log: 
        "results/extended/logs/hmmscan.txt"
    conda: 
        "hmmscan.yaml"
    resources: 
        load = 50
    params: 
        n=config['resources']["max_cpus"]
    shell: 
        "hmmsearch --cpu {params.n} \
        --domtblout {output} \
        results/extended/tr_2_prot/Pfam-A.hmm \
        {input.longest_orfs} \
        &> {log}"


rule transdecoder_Predict:
    """
    Predict coding status of longest ORFs
    """
    input: 
        transcripts_fa="results/extended/tr_2_prot/denovo.transcripts.fasta",
        pfam="results/extended/tr_2_prot/hmmscan/pfam.domtblout",
        pep="results/extended/tr_2_prot/longest_orfs.pep"
    output: 
        transdecoder_pep="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder.pep",
        transdecoder_gff3="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder.gff3",
        bed="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder.bed",
        cds="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder.cds"
    benchmark: 
        "results/extended/benchmarks/transdecoder_Predict.json"
    log: 
        "results/extended/logs/transdecoder_Predict.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config['resources']["max_cpus"],
        out_pep="denovo.transcripts.fasta.transdecoder.pep",
        out_gff3="denovo.transcripts.fasta.transdecoder.gff3",
        out_bed="denovo.transcripts.fasta.transdecoder.bed",
        out_cds="denovo.transcripts.fasta.transdecoder.cds"
    shell:
        "TransDecoder.Predict \
        -t {input.transcripts_fa} \
        --retain_pfam_hits {input.pfam} \
        --single_best_only \
        -O results/extended/tr_2_prot/ 2> {log} 1>&2 && \
        cp -f {params.out_pep} {output.transdecoder_pep} && \
        cp -f {params.out_gff3} {output.transdecoder_gff3} && \
        cp -f {params.out_bed} {output.bed} && \
        cp -f {params.out_cds} {output.cds} && \
        rm {params.out_pep} {params.out_gff3} {params.out_bed} {params.out_cds}"


rule gtf_to_alignment_gff3:
    input: 
        gff_denovo = "results/extended/tr_2_prot/denovo.transcripts.gtf",
    output: 
        denovo_alignment_gff3="results/extended/tr_2_prot/denovo_alignment.gff"
    benchmark: 
        "results/extended/benchmarks/gtf_to_alignment_gff3.json"
    log: 
        "results/extended/logs/gtf_to_alignment_gff3.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config['resources']["max_cpus"]
    shell: 
        "gtf_to_alignment_gff3.pl \
        {input.gff_denovo} > {output.denovo_alignment_gff3} 2> {log}"


rule transcript_sequences_alignment_orf_to_genome_orf:
    """
    Check if ORFs can be aligned to the genome
    """
    input: 
        transdecoder_gff3="results/extended/tr_2_prot/denovo.transcripts.fasta.transdecoder.gff3",
        denovo_alignment_gff3="results/extended/tr_2_prot/denovo_alignment.gff",
        denovo_fasta = "results/extended/tr_2_prot/denovo.transcripts.fasta"
    output: 
        "results/extended/tr_2_prot/transcripts.genome.gff3"
    benchmark: 
        "results/extended/benchmarks/transcript_sequences_alignment_orf_to_genome_orf.json"
    log: 
        "results/extended/logs/transcript_sequences_alignment_orf_to_genome_orf.txt"
    conda: 
        "transdecoder.yaml"
    params: 
        n=config['resources']["max_cpus"]
    shell: 
        "cdna_alignment_orf_to_genome_orf.pl {input.transdecoder_gff3} {input.denovo_alignment_gff3} {input.denovo_fasta} > {output} 2> {log}"


rule gff3_file_to_proteins:
    """
    Translate ORFs aligned to the genome
    """
    input: 
        transcripts_genome="results/extended/tr_2_prot/transcripts.genome.gff3",
        genome_fasta=branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/GRCh38.primary_assembly.genome_copy.fa',
                              otherwise = config["reference"]["genome_fasta"])
    output: 
        "results/extended/tr_2_prot/proteome_denovo.fasta"
    benchmark: 
        "results/extended/benchmarks/gff3_file_to_proteins.json"
    log: 
        "results/extended/logs/gff3_file_to_proteins.txt"
    conda: 
        "gff3_file_to_proteins.yaml"
    shell: 
        "gff3_file_to_proteins.pl \
        --gff3 {input.transcripts_genome} \
        --fasta {input.genome_fasta} > {output} \
        2> {log}"


rule rmdup_orf:
    """
    remove duplivate entries in the predicted orfs
    """
    input:
        proteome_denovo="results/extended/tr_2_prot/proteome_denovo.fasta"
    output:
        proteome_denovo_unique="results/extended/tr_2_prot/proteome_denovo_unique.fasta"
    benchmark: 
        "results/extended/benchmarks/rmdup_orf.json"
    log: 
        "results/extended/logs/rmdup_orf.txt"
    conda: 
        "R_remove_duplicates.yaml"
    params: 
        n=config['resources']["max_cpus"]
    script: 
        "remove_duplicate_proteins_ORF.R"


rule rmdup_denovo_and_reference_proteome:
    """
    keep only novel sequences, matched to the reference transcripts
    the others will be included as they are in a reference proteome 
    """
    input:
        rmdup_orf = "results/extended/tr_2_prot/proteome_denovo_unique.fasta",
        protein_all = config["reference"]["gencode_translation_sequences_fasta"]
    output:
        proteome_ref_denovo="results/extended/tr_2_prot/proteome_ref_denovo.fasta"
    benchmark: 
        "results/extended/benchmarks/rmdup_denovo_and_reference_proteome.json"
    log: 
        "results/extended/logs/rmdup_denovo_and_reference_proteome.txt"
    conda: 
        "R_remove_duplicates.yaml"
    params: 
        n=config['resources']["max_cpus"]
    script: 
        "remove_duplicate_proteins_denovo_and_reference_proteome.R"


rule tx2gene:
    """
    Create a gene-transcript mapping with reference and novel transcripts
    """
    input:
        gffcomp_tracking = "results/extended/stringtie_assemble/gffcompare/gffcmp.tracking",
        transcripts_fa = "results/extended/tr_2_prot/merged_reference.transcripts.fasta",
        gencode_annotation_gtf = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gtf',
                              otherwise = config["reference"]["reference_gtf"])
    output:
        tx2gene = "results/extended/tr_2_prot/tx2gene.csv"
    benchmark: 
        "results/extended/benchmarks/tx2gene.json"
    log: 
        "results/extended/logs/tx2gene.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    params: 
        n=config['resources']["max_cpus"]
    script: 
        "tx2gene.R"


rule expression_cutoffs_for_proteome_db_extended:
    """
    filters the proteome.denovo.fasta to keep the proteins with a minimum
    supporting number of reads in min_filt_samples
    """
    input:
        prot = "results/extended/tr_2_prot/proteome_ref_denovo.fasta",
        reference_gff3 = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gff3',
                              otherwise = config["reference"]["reference_gff3"]),
#        pred_prot = "results/RNAsamba_coding_potential/pred_prot.fasta",
        tx2gene = "results/extended/tr_2_prot/tx2gene.csv",
        gencode_metadata_TrEMBL = config["reference"]["gencode_metadata_TrEMBL"],
        gencode_metadata_SwissProt = config["reference"]["gencode_metadata_SwissProt"],
        tr_quant=expand("results/extended/trascriptome_quant/salmon/{sample}/quant.sf", sample = SAMPLES)
        #tr_quant=expand("results/trascriptome_quant/salmon_BAM/{sample}/quant.sf", sample = SAMPLES)
    output:
        expr_prot="results/extended/tr_2_prot/proteome_ref_denovo_expressed.fasta",
        gene_tr_prot="results/extended/tr_2_prot/gene_tr_prot.csv",
        gene_tr_prot_SP="results/extended/tr_2_prot/gene_tr_prot_SP.csv",
        transcript_expression="results/extended/trascriptome_quant/transcript_expression.csv",
        transcript_expression_tpm="results/extended/trascriptome_quant/transcript_expression_tpm.csv",
        gene_expression="results/extended/trascriptome_quant/gene_expression.csv",
        gene_expression_tpm="results/extended/trascriptome_quant/gene_expression_tpm.csv"
    benchmark: 
        "results/extended/benchmarks/expression_cutoffs_for_proteome_db.json"
    log: 
        "results/extended/logs/expression_cutoffs_for_proteome_db.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    resources: 
        load = 100
    params: 
        n=config['resources']["max_cpus"],
        min_filt_samples=config['program_specific']["min_filt_samples"],
        min_filt_counts=config['program_specific']["min_filt_counts"]
    script: 
        "expression_cutoffs_extended.R"


rule rmdup_RNAsamba:
    """
    removes the sequence duplicates from the final database
    """
    input:
        expr_prot_incl_deep="results/extended/tr_2_prot/proteome_ref_denovo_RNAsamba_expressed.fasta"
    output:
        expr_prot_incl_deep_nodup="results/extended/tr_2_prot/proteome_ref_denovo_RNAsamba_expressed_nodup.fasta"
    benchmark: 
        "results/extended/benchmarks/rmdup_RNAsamba.json"
    log: 
        "results/extended/logs/rmdup_RNAsamba.txt"
    conda: 
        "R_remove_duplicates.yaml"
    params: 
        n=config['resources']["max_cpus"]
    script: 
        "remove_duplicate_proteins_RNAsamba.R"


rule include_nORFs:
    input: 
        nORFs="data/reference/Mus_Musculus_nORFs.fasta",
        expr_prot="results/extended/tr_2_prot/proteome_ref_denovo_expressed.fasta",
        expr_prot_incl_deep_nodup="results/extended/tr_2_prot/proteome_ref_denovo_RNAsamba_expressed_nodup.fasta"
    output: 
        expr_prot_nORFs="results/extended/tr_2_prot/proteome_ref_denovo_nORFs_expressed.fasta",
        expr_prot_nORFs_incl_deep_nodup="results/extended/tr_2_prot/proteome_ref_denovo_RNAsamba_nORFs_expressed_nodup.fasta"
    benchmark: 
        "results/extended/benchmarks/include_nORFs.json"
    log: 
        "results/extended/logs/include_nORFs.txt"
    shell: 
        "cat {input.expr_prot} {input.nORFs} > {output.expr_prot_nORFs} && \
        cat {input.expr_prot_incl_deep_nodup} {input.nORFs} > {output.expr_prot_nORFs_incl_deep_nodup} \
        2> {log}"
