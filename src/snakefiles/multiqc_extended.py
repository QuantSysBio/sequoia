rule multiqc_extended:
    input: 
        gene_expression="results/extended/trascriptome_quant/gene_expression.csv",
        proteome_ref_denovo_expressed="results/extended/tr_2_prot/proteome_ref_denovo_expressed.fasta"
    output: 
        multiqc_report="results/extended/multiqc/multiqc_report.html",
    benchmark: 
        "results/extended/benchmarks/multiqc.txt"
    log: 
        "results/extended/logs/multiqc.txt"
    conda: 
        "multiqc.yaml"
    shell: 
        "multiqc -f results/extended/ \
        -o results/extended/multiqc \
         2> {log}"

rule wf_extended_done:
    input:
        "results/extended/multiqc/multiqc_report.html"
    output:
        workflow_out = "results/extended.txt"
    shell:
        'touch {output}'
