rule expression_cutoffs_for_proteome_db:
    """
    filters the proteome.denovo.fasta to keep the proteins with a minimum
    supporting number of reads in min_filt_samples
    """
    input:
        prot = features["reference"]["gencode_translation_sequences_fasta"],
        reference_gff3 = features["reference"]["reference_gff3"],
        gencode_metadata_TrEMBL = features["reference"]["gencode_metadata_TrEMBL"],
        gencode_metadata_SwissProt = features["reference"]["gencode_metadata_SwissProt"],
        tr_quant = expand("results/trascriptome_quant/salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        expr_prot = "results/tr_2_prot/proteome_expressed.fasta",
        gene_tr_prot = "results/tr_2_prot/gene_tr_prot.csv",
        gene_tr_prot_SP = "results/tr_2_prot/gene_tr_prot_SP.csv",
        transcript_expression = "results/trascriptome_quant/transcript_expression.csv",
        transcript_expression_tpm = "results/trascriptome_quant/transcript_expression_tpm.csv",
        gene_expression = "results/trascriptome_quant/gene_expression.csv",
        gene_expression_tpm = "results/trascriptome_quant/gene_expression_tpm.csv"
    benchmark: 
        "results/benchmarks/expression_cutoffs_for_proteome_db.json"
    log: 
        "results/logs/expression_cutoffs_for_proteome_db.txt"
    conda: 
        "expression_cutoffs_for_proteome_db.yaml"
    resources:
        ncpus = config["salmon_cpus"],
        mem = config["max_mem"],
        time = config["max_time"],
        load = 50
    params: 
        n = config["salmon_cpus"],
        min_filt_samples = config["min_filt_samples"],
        min_filt_counts = config["min_filt_counts"]
    script: 
        "expression_cutoffs.R"
