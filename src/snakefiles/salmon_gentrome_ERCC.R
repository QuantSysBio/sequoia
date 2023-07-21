### ---------------------------------------------- pre-process ERCC  ----------------------------------------------
# description:  Concatenate genome, transcriptome. Remove ERCC sequences from genome to avoid being treated as decoys
#               
# input:        - genome, transcriptome .fasta
#
# output:       
#               - gentrome.fasta
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)

input <- list()
input$ref_tr <- readDNAStringSet(snakemake@input[["transcriptome_fasta"]])
input$ref_genome <- readDNAStringSet(snakemake@input[["genome_fasta"]])

gentrome <- unique(c(input$ref_tr, input$ref_genome))
writeXStringSet(gentrome, unlist(snakemake@output[["gentrome"]]))
