rule multiqc:
    input:
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
        mem = config["resources"]["max_mem"]      
    shell: 
        """
        multiqc -f results/
        -o results/multiqc
        2> {log}
        """
