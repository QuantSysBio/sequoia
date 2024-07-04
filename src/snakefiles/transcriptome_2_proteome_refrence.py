rule expression_cutoffs_for_proteome_db:
    """
    filters the proteome.denovo.fasta to keep the proteins with a minimum
    supporting number of reads in min_filt_samples
    """
    input:
        prot = config["reference"]["gencode_translation_sequences_fasta"],
        reference_gff3 = branch(config["reference"]['ERCC_normalization'],
                              then = 'data/references/gencode.v45.primary_assembly.annotation_copy.gff3',
                              otherwise = config["reference"]["reference_gff3"]),
        gencode_metadata_TrEMBL = config["reference"]["gencode_metadata_TrEMBL"],
        gencode_metadata_SwissProt = config["reference"]["gencode_metadata_SwissProt"],
        tr_quant = expand("results/refrence/transcriptome_quant/salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        expr_prot = "results/refrence/tr_2_prot/proteome_expressed.fasta",
        gene_tr_prot = "results/refrence/tr_2_prot/gene_tr_prot.csv",
        gene_tr_prot_SP = "results/refrence/tr_2_prot/gene_tr_prot_SP.csv",
        transcript_expression = "results/refrence/transcriptome_quant/transcript_expression.csv",
        transcript_expression_tpm = "results/refrence/transcriptome_quant/transcript_expression_tpm.csv",
        gene_expression = "results/refrence/transcriptome_quant/gene_expression.csv",
        gene_expression_tpm = "results/refrence/transcriptome_quant/gene_expression_tpm.csv"
    benchmark: 
        "results/refrence/benchmarks/expression_cutoffs_for_proteome_db.json"
    log: 
        "results/refrence/logs/expression_cutoffs_for_proteome_db.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    resources:
        ncpus = config["program_specific"]["salmon_cpus"],
        mem = config["resources"]["max_mem"],
        load = 50
    params: 
        n = config["program_specific"]["salmon_cpus"],
        min_filt_samples = config["program_specific"]["min_filt_samples"],
        min_filt_counts = config["program_specific"]["min_filt_counts"]
    script: 
        "expression_cutoffs.R"
