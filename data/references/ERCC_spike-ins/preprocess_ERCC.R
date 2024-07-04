### ---------------------------------------------- pre-process ERCC  ----------------------------------------------
# description:  Modify ERCC fasta and gtf attributes to fit GENCODE header system 
#               
# input:        1. ERCC .fasta
#               2. ERCC .gtf
# output:       
#               - fasta, gtf and gff files.
#               - Concatentate them with genome, transcriptome and transcriptome annotation files
#               
# author:       YH

library(Biostrings)
library(plyranges)
library(dplyr)

input <- list()
input$gtf <- read_gff("data/reference/ERCC_spike-ins/ercc.gtf")
input$fasta <- readDNAStringSet("data/reference/ERCC_spike-ins/ercc.fa")

# Add transcript_id attribute to ERCC GTF
gtf <- input$gtf %>%
  mutate(input$gtf, transcript_id = seqnames) %>%
  as_tibble() %>%
  select(c(colnames(as_tibble(input$gtf)), "transcript_id")) %>%
  as(Class = "GRanges")

gtf <- c(gtf, input$gtf) %>%
  arrange(seqnames)

fasta <- input$fasta
names(fasta)[1]

# Modify fasta header accordingly
fasta_names <- names(fasta) %>%
  str_split_fixed(pattern = " ", n = 4) %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(transcript_id = V1) %>%
  left_join(as_tibble(select(gtf, transcript_id, gene_id))) %>%
  mutate(header = paste(transcript_id, gene_id, V2, V3, V4, sep="|"))
fasta_names

names(fasta) <- fasta_names$header

# Save outputs
write_gff2(gtf, "data/reference/ERCC_spike-ins/ercc_gencode.gtf")
write_gff3(gtf, "data/reference/ERCC_spike-ins/ercc_gencode.gff")
writeXStringSet(fasta, "data/reference/ERCC_spike-ins/ercc_gencode.fasta")


# Add ERCC to reference and export. (change filenames if necessary)
input$ref_gtf <- read_gff("data/reference/gencode.v38.primary_assembly.annotation.gtf.gz")
input$ref_gff <- read_gff("data/reference/gencode.v38.primary_assembly.annotation.gff3.gz")
input$ref_tr <- readDNAStringSet("data/reference/gencode.v38.transcripts.fa.gz")
input$ref_genome <- readDNAStringSet("data/reference/GRCh38.primary_assembly.genome.fa.gz")

ref_gtf <- c(input$ref_gtf, gtf)
ref_gtf
ref_tr <- c(input$ref_tr, fasta)
ref_tr

# Concatenation order matters for Salmon. Make sure that in gentrome .fasta ERCC are preceeding genomic sequences
ref_genome <- c(fasta, input$ref_genome)
ref_genome

write_gff2(ref_gtf, "data/reference/gencode.v38.primary_assembly.annotation.gtf")
write_gff3(ref_gtf, "data/reference/gencode.v38.primary_assembly.annotation.gff3")
writeXStringSet(ref_tr, "data/reference/gencode.v38.transcripts.fa")
writeXStringSet(ref_genome, "data/reference/GRCh38.primary_assembly.genome.fa")

print("DONE")

