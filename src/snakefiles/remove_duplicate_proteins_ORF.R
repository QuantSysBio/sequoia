### ---------------------------------------------- Remove duplicates  ----------------------------------------------
# description:  Remove duplicates
#               
# input:        1. De-novo proteome
# output:       
#               - Unique de-novo proteome, duplicate sequences are annotated as in reference
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)
input <- readAAStringSet(snakemake@input[["proteome_denovo"]])

out <- unique(input)

writeXStringSet(out, unlist(snakemake@output[["proteome_denovo_unique"]]))
