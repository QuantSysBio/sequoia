### SPI-art project ###
# 
# description:  In silico translate transcriptome CDS, match to Gencode translation sequences. 6-frame translation and tags
# input:        - human Gencode transcriptome annotation .gtf
#               - human Gencode genome package 
# output:       A folder with: .gtf file of selected features, .fasta of 6-frame translation.
#               132 sequences don't match the genetic code and will be replaced by Gencode sequences
# author:       Yehor Horokhovskyi




library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(plyranges)
library(stringr)
library(readr)

###-------------------Creating the result folder--------------------------
create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
}
folder_path <- 'results/Exhaustive/'
create_folder_if_not_exists(folder_path)
### ------------------ Start user variables ------------------------------

# Minimal length filter
min_ORF_length = snakemake@params[['min_orf_length']]

gff_path = snakemake@input[['reference_gff3']]
transcripts_path = snakemake@input[['reference_transcripts']]
prot_path = snakemake@input[['reference_prot']]
gff <- read_gff(gff_path)
txdb <- makeTxDbFromGFF(file = paste0(gff_path),
                        dataSource="GENCODE_v45",
                        organism="Homo sapiens")
#txdb

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
# Hsapiens <- readRDS("data/BSgenome.Hsapiens.UCSC.hg38.rds")

# update the names to match
upd_names <- names(Hsapiens) %>% 
  str_replace_all(pattern = "v1", replacement = "\\.1") %>% 
  str_replace_all(pattern = "v2", replacement = "\\.2") %>% 
  str_extract(pattern = str_c(unique(as.character(seqnames(gff)))[26:length(unique(as.character(seqnames(gff))))], collapse = "|"))
seqnames(Hsapiens)[which(!is.na(upd_names))] <- upd_names[which(!is.na(upd_names))]

# Load expression data
# expr <- vroom(paste0("data/expressed_all.csv"))

# Load GENCODE
# Substitute selencysteine with cysteine and remove X
Gencode_prot <- readAAStringSet(paste0(prot_path)) %>%
  # stringr::str_replace_all(pattern = "[[X]]", replacement = "") %>%
  # stringr::str_replace_all(pattern = "U", replacement = "C") %>%
  AAStringSet()
# names(Gencode_prot) <- names(readAAStringSet(paste0(dir_K562, "data/gencode.v33.pc_translations.fa")))
Gencode_prot

Gencode_RNA <- readDNAStringSet(paste0(transcripts_path)) 
Gencode_RNA

### ------------------ End user variables ------------------------------
# Match Gencode mRNA with txDB transcript sequences
CDS_seq <- extractTranscriptSeqs(Hsapiens, txdb)
table(Gencode_RNA %in% CDS_seq)

# Extract CDS sequences - they don't match the full mRNAs
CDS <- cdsBy(txdb, use.names=TRUE)
CDS_seq <- extractTranscriptSeqs(Hsapiens, CDS)
table(Gencode_RNA %in% CDS_seq)

# ---------------------------  Translate --------------------------- 
# Demo: 6-frame translation matches almost all sequences to Gencode reference proteins
CDS_subseqs <- lapply(1:3, function(pos) 
  subseq(c(DNAStringSet(CDS_seq), reverseComplement(DNAStringSet(CDS_seq))), start=pos))

# translate 6-frame
CDS_seq_6_frame <- lapply(CDS_subseqs, translate, if.fuzzy.codon = c("solve"), no.init.codon = T)
CDS_seq_6_frame <- c(CDS_seq_6_frame[[1]], 
                     CDS_seq_6_frame[[2]], 
                     CDS_seq_6_frame[[3]])
names(CDS_seq_6_frame) <- paste0(rep(c("forward_f_1|", "complement_f_1|", 
                                       "forward_f_2|", "complement_f_2|", 
                                       "forward_f_3|", "complement_f_3|"), 
                                     each = length(CDS_seq)), 
                                 names(CDS_seq_6_frame))
CDS_seq_6_frame


# NT names
CDS_subseq_names <- paste0(rep(c("forward_f_1|", "complement_f_1|", 
                                 "forward_f_2|", "complement_f_2|", 
                                 "forward_f_3|", "complement_f_3|"), 
                               each = length(CDS_subseqs[[1]]) / 2), 
                           names(CDS_subseqs[[1]]))

CDS_subseqs <- unlist(DNAStringSetList(CDS_subseqs))
names(CDS_subseqs) <- CDS_subseq_names

# Compare CDS & off-frame
tmp <- names(CDS_seq_6_frame)
CDS_seq_6_frame <- gsub(x = CDS_seq_6_frame, pattern = "U", replacement = "C") %>%
  # stringr::str_replace_all(pattern = "[[*]]", replacement = "") %>%
  # stringr::str_replace_all(pattern = "[[X]]", replacement = "") %>%
  # gsub(pattern = "U", replacement = "C") %>%
  AAStringSet()
names(CDS_seq_6_frame) <- tmp

CDS_seq_6_frame
x <- Gencode_prot %in% CDS_seq_6_frame
table(x)
Gencode_prot[x]
Gencode_prot[!Gencode_prot %in% CDS_seq_6_frame]
x <- CDS_seq_6_frame %in% Gencode_prot
table(x)
CDS_seq_6_frame[x]
CDS_seq_6_frame[!CDS_seq_6_frame %in% Gencode_prot]

# ---------------------------  Match to Gencode --------------------------- 
### Create a 6-frame translation DB where the reference sequences come from Gencode + all other frames
# Reference sequences inherit the Gencode tags
# Replace the matching sequences with Gencode

CDS_annot <- AAStringSet()
CDS_annot_NT <- DNAStringSet()
seq_along(Gencode_prot)[1:2]
for (i in seq_along(Gencode_prot)[1:2]) {
  # i = 1
  #rm(char_dist)
  print(i)
  print('c-----------------------------------------------')
  # Select translation sequences of the same transcript
  x <- Gencode_prot[i]
  print(x)
  xn <- x %>%
    names() %>%
    str_split_fixed("[[|]]", Inf) %>%
    as.data.frame()
  print(xn)
  
  y <- CDS_seq_6_frame[grep(xn$V2, names(CDS_seq_6_frame))]
  print(y)
  print("Done")
  
  # If the exact matching doesn't work, match to the closest 6-frame translation from the same transcript
  if (x %in% y) {
    
    char_dist <- adist(x, y)
    yn <- names(y)[which.min(char_dist)] %>%
      paste0("|", "main_ORF","|mismatches:", min(char_dist), "|",names(x))
    
    yn_2 <- names(y)[-which.min(char_dist)] %>%
      paste0("|", "off-frame","|mismatches:", char_dist[-which.min(char_dist)])
    
    # Create output AAStringSet
    out <- AAStringSet(x = x)
    names(out) <- yn
    
    out_2 <- AAStringSet(x = y[-which.min(char_dist)])
    names(out_2) <- yn_2
    out <- c(out, out_2)
  } else {
    
    # Partial match
    char_dist <- adist(x, y)
    yn <- names(y)[which.min(char_dist)] %>%
      paste0("|", "main_ORF","|mismatches:", min(char_dist), "|",names(x))
    
    yn_2 <- names(y)[which(!char_dist %in% which.min(char_dist))] %>%
      paste0("|", "off-frame","|mismatches:", char_dist[which(!char_dist %in% which.min(char_dist))])
    
    # Create output AAStringSet
    out <- AAStringSet(x = x)
    names(out) <- yn
    
    out_2 <- AAStringSet(x = y[which(!char_dist %in% which.min(char_dist))])
    names(out_2) <- yn_2
    
    out <- c(out, out_2)
  }
  CDS_annot <- c(CDS_annot, out)
}
# CDS_annot %>% names() %>% View()
print(min_ORF_length)
# ---------------------------  Filter by expression & length --------------------------- 
print(width(CDS_annot))
print(rep(min_ORF_length, length(CDS_annot)))
print(width(CDS_annot) >= rep(min_ORF_length, length(CDS_annot)))
print(CDS_annot[width(CDS_annot) >= rep(min_ORF_length, length(CDS_annot))])
CDS_annot <- CDS_annot[width(CDS_annot) >= rep(min_ORF_length, length(CDS_annot))]

CDS_annot_names <- CDS_annot %>% 
  names() %>%
  str_split_fixed("[[|]]", Inf) %>%
  as.data.frame()

# CDS_expressed <- CDS_annot[CDS_annot_names$V2 %in% expr$TXNAME[expr$tr_expressed == TRUE]]

# ---------------------------  Save fasta --------------------------- 
CDS_annot[str_detect(names(CDS_annot), "main_ORF")] %>%
  unique() %>%
  writeXStringSet(paste0(snakemake@output[['CDS_main_ORF_fasta']]))

CDS_annot[str_detect(names(CDS_annot), "off-frame")] %>%
  unique() %>%
  writeXStringSet(paste0(snakemake@output[['CDS_frameshift_fasta']]))

CDS_subseqs %>%
  writeXStringSet(paste0(snakemake@output[['CDS_NT_fasta']]))



