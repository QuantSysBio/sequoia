### SPI-art project ###
# 
# description:  In silico translate transcriptome UTRs, filter by expression, translate in 3 frames
# input:        - human Gencode transcriptome annotation .gtf
#               - human Gencode genome package 
# output:       A folder with: .fasta of 3-frame translation.
#               
# author:       Yehor Horokhovskyi

library(bettermc)
library(BSgenome.Hsapiens.GRCh38.primary.assembly)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(tidyr)
library(ORFik)  
library(plyranges)
library(stringr)
library(vroom)

### ------------------ Start user variables ------------------------------
dir.create(paste0("results"))
dir.create(paste0("results/all"))
dir.create(paste0("results/expressed"))
dir.create(paste0("results/curated"))

# Multi-threading
Ncpu = 7

# Genome object
Hsapiens <- BSgenome.Hsapiens.GRCh38.primary.assembly
Hsapiens

# For non-CDS ORF discovery
min_ORF_length = 8

# Load expression data
expr <- vroom(paste0("data/expressed_all.csv"))
expr_tr <- expr %>%
  filter(tr_expressed == TRUE) %>%
  pull(TXNAME) %>%
  unique()

# Load nORFs
nORFs <- read_gff("data/nORFs/nORFsDB.1.1.gtf")
nORFs_annot <- vroom("data/nORFs/nORFsDB.1.1.classification.tsv")

# Load OpenProt
OpenProt <- readAAStringSet("data/OpenProt/human-openprot-r1_6-altprots+isoforms_min_1_pep-grch38.95.fasta")

# Load SmProt
SmProt <- list.files("data/smProt/", pattern = ".fa.gz", full.names = T) %>%
  lapply(readAAStringSet) %>%
  AAStringSetList() %>%
  unlist()

### ------------------ End user variables ------------------------------
# nORFs
nORFs_all <- nORFs %>%
  as_tibble() %>%
  cbind(nORFs_annot) %>%
  as_tibble()

nORFs_fasta <- AAStringSet(nORFs_all$AA_seq)
names(nORFs_fasta) <- nORFs_all$novelORF_ID
nORFs_fasta <- unique(nORFs_fasta)
nORFs_fasta <- nORFs_fasta[width(nORFs_fasta) >= min_ORF_length]

writeXStringSet(nORFs_fasta, "results/curated/nORFs.fasta", append = F)

# OpenProt
OpenProt <- OpenProt[width(OpenProt) >= min_ORF_length]
OpenProt <- unique(OpenProt)
writeXStringSet(OpenProt, "results/curated/OpenProt.fasta", append = F)

# SmProt
SmProt <- SmProt[width(SmProt) >= min_ORF_length]
SmProt <- unique(SmProt)
writeXStringSet(SmProt, "results/curated/SmProt.fasta", append = F)
