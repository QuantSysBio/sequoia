### SPI-art project ###
# 
# description:  Create a BSgenome for a Gencode genome used
# input:        - human Gencode genome .fasta.gz
#               - primary_assembly-seed (see BSgenome package)
# 
# output:       An R package BSgenome.Hsapiens.GRCh38.primary.assembly_1.0.tar.gz
# 
# author:       Yehor Horokhovskyi

library(Biostrings)
library(BSgenome)
library(stringr)

### This file contains information about genome package
seed_files <- list.files("data/",pattern = "-seed",
                         full.names = T)
seqs_srcdir <- readLines(seed_files)
cat(seqs_srcdir, sep="\n")
seqs_srcdir <- seqs_srcdir[str_starts(seqs_srcdir, "seqs_srcdir:")]
seqs_srcdir <- str_split_fixed(seqs_srcdir, "seqs_srcdir: ", 2)[,2]
suppressWarnings(dir.create(seqs_srcdir))

### Split the fa.gz gencode genome 
genome <- readDNAStringSet("data/GRCh38.primary_assembly.genome.fa.gz")
for (i in seq_along(genome)) {
  chromosome <- genome[i]
  chromosome_name <- str_split_fixed(names(genome)[[i]], " ", 2)[,1]
  writeXStringSet(chromosome, filepath = paste0(seqs_srcdir, chromosome_name, ".fa"))
  rm(chromosome, chromosome_name)
}

###  Create and install a package
# destdir <- str_replace(seqs_srcdir, ".fa.gz.split", ".package")
destdir <- str_remove(seqs_srcdir, ".genome.fa.gz.split")

unlink(destdir, recursive = T, force = T)
dir.create(destdir)
forgeBSgenomeDataPkg(x = seed_files, 
                     verbose = TRUE,
                     seqs_srcdir = seqs_srcdir, 
                     destdir = destdir)

system(paste0("R CMD build data/GRCh38.primary_assembly/BSgenome.Hsapiens.GRCh38.primary.assembly"))
system(paste0("R CMD check BSgenome.Hsapiens.GRCh38.primary.assembly_1.0.tar.gz"))

install.packages("BSgenome.Hsapiens.GRCh38.primary.assembly_1.0.tar.gz", 
                 repos = NULL, type = "source")

library(BSgenome.Hsapiens.GRCh38.primary.assembly)
Hsapiens <- BSgenome.Hsapiens.GRCh38.primary.assembly
