### ---------------------------------------------- Remove redundant transcripts  ----------------------------------------------
# description:  Create a hybrid transcriptome annotation from reference and stringtie output.
#               Keep denovo transcripts that are not equivalent according to gffcompare. 
#               
# input:        1. Reference and de-novo transcriptomes (.gtf and .fasta)
# output:       
#               - de-novo only transcripts (.fasta and .gtf2)
#               - combined reference and de-novo transcripts (.fasta and .gff3)
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)
library(dplyr)
library(plyranges)

ref.fa <- readDNAStringSet(snakemake@input[["ref_fa"]])
mer.fa <- readDNAStringSet(snakemake@input[["mer_fa"]])

ref.gtf <- read_gff(snakemake@input[["ref_gtf"]])
mer.gtf <- read_gff(snakemake@input[["mer_gtf"]])

# Select novel transcripts
new.tr.id <- mer.gtf$transcript_id[!mer.gtf$class_code == "="] %>% 
  na.omit() %>%
  unique()

new.gtf <- mer.gtf[mer.gtf$transcript_id %in% new.tr.id,]
new.fa <- mer.fa[names(mer.fa) %in% new.tr.id]

# Combine with reference
combined.gtf <- c(ref.gtf, new.gtf)
combined.fa <- c(ref.fa, new.fa)

# Export
writeXStringSet(combined.fa, unlist(snakemake@output[["fasta"]]))
writeXStringSet(new.fa, unlist(snakemake@output[["denovo_fasta"]]))

write_gff3(combined.gtf, unlist(snakemake@output[["gff_final"]]))
write_gff2(new.gtf, unlist(snakemake@output[["gtf_denovo"]]))
