rule multiqc:
    input:
        #proteome_ref_denovo_expressed="results/tr_2_prot/proteome_ref_denovo_expressed.fasta",
        gene_expression="results/trascriptome_quant/gene_expression.csv"
    output: 
        multiqc_report="results/multiqc/multiqc_report.html"
    benchmark: 
        "results/benchmarks/multiqc.txt"
    log: 
        "results/logs/multiqc.txt"
    conda: 
        "multiqc.yaml"
    resources:
        ncpus = 1,
        mem = config["max_mem"],
        time = config["max_time"]        
    shell: 
        "multiqc -f results/ \
        -o results/multiqc \
         2> {log}"

