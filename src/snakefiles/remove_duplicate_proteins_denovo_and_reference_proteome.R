### ---------------------------------------------- Remove duplicates  ----------------------------------------------
# description:  Remove duplicates
#               
# input:        1. Reference proteome
#               2. De-novo proteome
# output:       
#               - Unique proteome, duplicate sequences are annotated as in reference
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)
rmdup_orf <- readAAStringSet(snakemake@input[["rmdup_orf"]])
protein_all <- readAAStringSet(snakemake@input[["protein_all"]])

out <- unique(c(protein_all, rmdup_orf))

writeXStringSet(out, unlist(snakemake@output[["proteome_ref_denovo"]]))
