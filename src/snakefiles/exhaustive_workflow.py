rule CDS:
    """
    In silico translate transcriptome CDS, match to Gencode translation sequences. 6-frame translation and tags
    """
    input:
        reference_gff3 = config["reference"]["reference_gtf_no_ERCC"],
        reference_transcripts = config["reference"]["transcriptome_fasta_no_ERCC"],
        reference_prot = config['reference']['gencode_translation_sequences_fasta_no_ERCC']
    output:
        CDS_main_ORF_fasta = "results/Exhaustive/CDS/CDS_main_ORF.fasta",
        CDS_frameshift_fasta = "results/Exhaustive/CDS/CDS_frameshift.fasta",
        CDS_NT_fasta = "results/Exhaustive/CDS/CDS_NT.fasta"
    benchmark: 
        "results/Exhaustive/benchmarks/CDS.json"
    log: 
        "results/Exhaustive/logs/CDS.txt"
    conda: 
        "01_CDS.yaml"
    params: 
        min_orf_length = config["program_specific"]["min_orf_length"]
    script: 
        "01_CDS.R"


rule cryptic:
    """
    In silico translate transcriptome CDS, match to Gencode translation sequences. 6-frame translation and tags
    """
    input:
        reference_gff3 = config["reference"]["reference_gtf_no_ERCC"],
        CDS_out = "results/Exhaustive/CDS/CDS_NT.fasta"
    output:
        all_result = directory("results/Exhaustive/cryptic/all"),
        Workflow_output = "results/exhaustive.txt"
    benchmark: 
        "results/Exhaustive/benchmarks/cryptic.json"
    log: 
        "results/Exhaustive/logs/cryptic.txt"
    conda: 
        "02_cryptic_strata.yaml"
    params:
        init_codons = config['program_specific']['start_codons']
    script: 
        "02_cryptic_strata.R"
