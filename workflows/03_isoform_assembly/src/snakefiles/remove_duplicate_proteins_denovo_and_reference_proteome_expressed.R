### ---------------------------------------------- Remove duplicates  ----------------------------------------------
# description:  Remove duplicates
#               
# input:        Expressed proteome
# output:       
#               - Unique expressed proteome
#               
# author:       YH


### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)
expr_prot <- readAAStringSet(snakemake@input[["expr_prot"]])

out <- unique(expr_prot)

writeXStringSet(out, unlist(snakemake@output[["expr_prot_nodup"]]))
