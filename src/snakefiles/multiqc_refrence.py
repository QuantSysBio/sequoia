rule multiqc:
    input:
        gene_expression="results/refrence/transcriptome_quant/gene_expression.csv"
    output: 
        multiqc_report="results/refrence/multiqc/multiqc_report.html"
        
    benchmark: 
        "results/refrence/benchmarks/multiqc.txt"
    log: 
        "results/refrence/logs/multiqc.txt"
    conda: 
        "multiqc.yaml"
    resources:
        ncpus = 1,
        mem = config["resources"]["max_mem"]      
    shell: 
        "multiqc -f results/refrence/ -o results/refrence/multiqc 2> {log}"

rule wf_reference_done:
    input:
        "results/refrence/multiqc/multiqc_report.html"
    output:
        workflow_out = "results/reference.txt"
    shell:
        'touch {output}'
