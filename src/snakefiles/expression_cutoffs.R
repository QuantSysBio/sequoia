### ---------------------------------------------- Expression cutoffs  ----------------------------------------------
# description:  Normalize transcript expression and create expressed proteome
#               
# input:        1. Proteomes of interest
#               2. Gene-transcript mappingg
#               3. Expression cutoffs
#               4. Salmon quantitation of RNA expression
# output:       
#               - Gene-transcript-protein mapping, transcript expression, expressed proteomes .fasta
#               
# author:       YH

library(DESeq2)
library(dplyr)
library(plyranges)
library(readr)
library(seqinr)
library(stringr)
library(tximport)

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = c("output", "message"))
  
  # load the proteome
  prot <- read.fasta(snakemake@input[["prot"]], whole.header=T)
  
  # transcript-level filtering options
  min_filt_samples <- snakemake@params[["min_filt_samples"]] # in at leat # samples
  min_filt_counts <- snakemake@params[["min_filt_counts"]]  # min counts 
  
  # Load the GENCODE mapping transcript-Uniprot
  tr_SwPr_colnames <- c("TXNAME", "SwissProt", "SwissProt_2")
  tr_SwPr <- read_tsv(snakemake@input[["gencode_metadata_SwissProt"]], col_names = tr_SwPr_colnames)
  tr_SwPr$SwissProt_2 <- NULL
  
  # Load the GENCODE mapping transcript-TrEMBL
  tr_TrEMBL_colnames <- c("TXNAME", "TrEMBL", "TrEMBL_2")
  tr_TrEMBL <- read_tsv(snakemake@input[["gencode_metadata_TrEMBL"]], col_names = tr_TrEMBL_colnames)
  tr_TrEMBL$TrEMBL_2 <- NULL
  
  # Salmon inputs
  reference_gff3 <- read_gff(snakemake@input[["reference_gff3"]])
  expr.files <- as.character(snakemake@input[["tr_quant"]])
}

# Salmon inputs
expr.names <- expr.files %>%
  str_remove_all(pattern = "\\/quant.sf") %>%
  str_split_fixed(pattern = "[/]", n = Inf) %>%
  as_tibble()

expr.names <- pull(expr.names, ncol(expr.names)) 
expr.names

tx2gene <- reference_gff3 %>%
  as_tibble() %>%
  filter(type == "transcript") %>%
  select(transcript_id, gene_id) %>%
  unique() %>%
  dplyr::rename(GENEID = gene_id, 
         TXNAME = transcript_id)


# --------------------------------------------------- End user inputs --------------------------
# Import abundance estimates on transcript level
import <- tximport(files = expr.files, 
                   importer = data.table::fread,
                   type = "salmon", 
                   tx2gene = tx2gene, 
                   countsFromAbundance = "dtuScaledTPM", 
                   ignoreTxVersion = F, 
                   ignoreAfterBar = T, 
                   txOut = T)

print('1111111111')
colnames(import$counts) <- expr.names
colnames(import$abundance) <- expr.names

### Consider ERCCs present if at least 50% of sequences pass the expression cutoffs
ERCC_present <- any(grepl("ERCC-", rownames(import[["counts"]])))
ERCC_counts <- import[["counts"]][grepl("ERCC-", rownames(import[["counts"]])),]
if (dim(import[["counts"]])[[2]] == 1) {
  ERCC_present <- sum((ERCC_counts >= min_filt_counts) >= min_filt_samples) / length(ERCC_counts) > 0.5
} else {
  ERCC_present <- sum(rowSums(ERCC_counts >= min_filt_counts) >= min_filt_samples) / nrow(ERCC_counts) > 0.5
}
if (is.na(ERCC_present)) {
  ERCC_present <- FALSE
}

if (ERCC_present) {
  # If ERCC sequences are present, use them as internal controls
  condition <- data.frame(sample = expr.names)
  
  if (length(expr.names) == 1) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import$counts), colData = condition, ~ 1)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import$counts), colData = condition, ~ sample)
  }
  dds <- DESeq2::estimateSizeFactors(dds, controlGenes=grep("ERCC-", rownames(import[["counts"]])))
  ERCC_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  # Select expressed transcripts
  if (dim(import[["counts"]])[[2]] == 1) {
  idx <- (ERCC_counts >= min_filt_counts ) >= min_filt_samples
} else {
  idx <- rowSums(ERCC_counts >= min_filt_counts ) >= min_filt_samples
}
#  idx <- (ERCC_counts >= min_filt_counts ) >= min_filt_samples
  expressed_tr <- ERCC_counts[idx,] %>% as.data.frame()
  
  tmp <- rownames(expressed_tr)
  tmp <- stringr::str_split_fixed(tmp, pattern="\\|",n=Inf) %>% as.data.frame()
  expressed_tr <- cbind(tmp$V1, expressed_tr)
  colnames(expressed_tr)[1] <- "TXNAME"
  
  # Normalize TPM using known lib.size factors
  expressed_tr_tpm <- sweep(as.data.frame(import$abundance[idx,]), 2, dds$sizeFactor, FUN = '/') %>% 
    as.data.frame()
  
  tmp <- rownames(expressed_tr)
  tmp <- stringr::str_split_fixed(tmp, pattern="\\|",n=Inf) %>% as.data.frame()
  expressed_tr_tpm <- cbind(tmp$V1, expressed_tr_tpm)
  colnames(expressed_tr_tpm)[1] <- "TXNAME"
  rownames(expressed_tr_tpm) <- 1:length(rownames(expressed_tr_tpm))
  
} else {
  condition <- data.frame(sample = expr.names)
  
  if (length(expr.names) == 1) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import$counts), colData = condition, ~ 1)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import$counts), colData = condition, ~ sample)
  }
  dds <- DESeq2::estimateSizeFactors(dds)

  #idx <- rowSums(import$counts >= min_filt_counts ) >= min_filt_samples
  if (dim(import[["counts"]])[[2]] == 1) {
  idx <- (import$counts >= min_filt_counts ) >= min_filt_samples
} else {
  idx <- rowSums(import$counts >= min_filt_counts ) >= min_filt_samples
}
  expressed_tr <- import$counts[idx,] %>% as.data.frame()
  
  tmp <- rownames(expressed_tr)
  tmp <- stringr::str_split_fixed(tmp, pattern="\\|",n=Inf) %>% as.data.frame()
  expressed_tr <- cbind(tmp$V1, expressed_tr)
  colnames(expressed_tr)[1] <- "TXNAME"
  
  # Export the transcript-level expression
  expressed_tr_tpm <- import$abundance[idx,] %>% as.data.frame()
  tmp <- rownames(expressed_tr)
  tmp <- stringr::str_split_fixed(tmp, pattern="\\|",n=Inf) %>% as.data.frame()
  expressed_tr_tpm <- cbind(tmp$V1, expressed_tr_tpm)
  colnames(expressed_tr_tpm)[1] <- "TXNAME"
  rownames(expressed_tr_tpm) <- 1:length(rownames(expressed_tr_tpm))
}

# ------------------- Create transcript to protein mapping
# 1. Ensembl proteins
Ensembl_proteins <- names(prot)
Ensembl_proteins <- names(prot)[grep(pattern = "ENS", x = Ensembl_proteins)]
Ensembl_proteins_names <- stringr::str_split_fixed(Ensembl_proteins, pattern="\\|",n=Inf) %>% as.data.frame()
EnsP <- data.frame("protein" = Ensembl_proteins_names$V1,
                   "TXNAME" = Ensembl_proteins_names$V2)

# 2. de-novo transdecoder proteins
Ensembl_proteins <- names(prot)
transdecoder_proteins <- names(prot)[-grep(pattern = "ENS", x = Ensembl_proteins)]
transdecoder_proteins_names <- stringr::str_split_fixed(transdecoder_proteins, pattern=" ",n=Inf) %>% as.data.frame()
transdecoder_proteins_names2 <- stringr::str_split_fixed(transdecoder_proteins_names$V1 , pattern=".p",n=Inf) %>% as.data.frame()
Tcons_P <- data.frame("protein"=transdecoder_proteins_names$V1,
                      "TXNAME"=transdecoder_proteins_names2$V1)

tx2prot <- rbind(EnsP, Tcons_P)

# ------------------- Keep the expressed proteins
expr_prot <- left_join(expressed_tr, tx2prot)
expr_prot <- left_join(expr_prot, tr_SwPr)
expr_prot <- left_join(expr_prot, tr_TrEMBL)


Ensembl_proteins <- names(prot)
Ensembl_proteins_names <- stringr::str_split_fixed(Ensembl_proteins, pattern="\\|",n=Inf) %>% as.data.frame()

which_ENSP_TCONS <- stringr::str_split_fixed(string = Ensembl_proteins_names$V1, pattern=" ",n=Inf) %>% as.data.frame()
which_ENSP_TCONS <- which_ENSP_TCONS$V1 %in% expr_prot$protein
expr_prot_fasta <- prot[which_ENSP_TCONS]

### Create a common mapping for all genes, transcripts and proteins
gene_tr_prot <- left_join(tx2gene, tx2prot) %>% unique()
gene_tr_prot <- data.frame("GENEID"=gene_tr_prot$GENEID,
                           "TXNAME"=gene_tr_prot$TXNAME,
                           "protein"=gene_tr_prot$protein)

# Add Swiss-Prot annotation from GENCODE
# idk where the unmapped transcripts are, because they are also not present in GENCODE transcriptome
gene_tr_prot_SP <- left_join(gene_tr_prot, tr_SwPr)
gene_tr_prot_SP$SwissProt_2 <- NULL
gene_tr_prot_SP <- left_join(gene_tr_prot_SP, tr_TrEMBL)

gene_tr_prot_SP <- data.frame("GENEID"=gene_tr_prot_SP$GENEID,
                              "TXNAME"=gene_tr_prot_SP$TXNAME,
                              "protein"=gene_tr_prot_SP$protein,
                              "SwissProt"=gene_tr_prot_SP$SwissProt,
                              "TrEMBL"=gene_tr_prot_SP$TrEMBL)

length(unique(tr_SwPr$SwissProt))
length(unique(gene_tr_prot_SP$SwissProt))
unmapped <- tr_SwPr[!tr_SwPr$SwissProt %in% gene_tr_prot_SP$SwissProt,]

expr_prot_tpm <- left_join(expressed_tr_tpm, tx2prot)
expr_prot_tpm <- left_join(expr_prot_tpm, tr_SwPr)
expr_prot_tpm <- left_join(expr_prot_tpm, tr_TrEMBL)

# --------------------------------------- Gene-level expression -----------------------------------------------------
rm(import)
gc()

# Import abundance estimates on gene level
import.gene <- tximport(files = expr.files, 
                        importer = data.table::fread,
                        type = "salmon", 
                        tx2gene = tx2gene, 
                        countsFromAbundance = "lengthScaledTPM", 
                        ignoreTxVersion = F, 
                        ignoreAfterBar = T, 
                        txOut = F)

print('22222')
if (ERCC_present) {
  
  # If ERCC sequences are present, use them as internal controls
  condition <- data.frame(sample = expr.names)
  
  
  if (length(expr.names) == 1) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import.gene$counts), colData = condition, ~ 1)
  } else {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(import.gene$counts), colData = condition, ~ sample)
  }
  dds <- DESeq2::estimateSizeFactors(dds, controlGenes=grep("ERCC-", rownames(import.gene[["counts"]])))
  
  ERCC_counts_gene <- DESeq2::counts(dds, normalized=TRUE)
  
  # Export the gene-level expression
  gene_expression <- data.frame("GENEID"=rownames(ERCC_counts_gene), 
                                ERCC_counts_gene %>% as.data.frame()) %>% as.data.frame()
  rownames(gene_expression) <- 1:length(rownames(gene_expression))
  colnames(gene_expression) <- c("GENEID", expr.names)
  
  
  gene_expression.tpm <- data.frame("GENEID"=rownames(import.gene$abundance), 
                                    sweep(as.data.frame(import.gene$abundance), 2, dds$sizeFactor, FUN = '/')) %>% as.data.frame()
  rownames(gene_expression.tpm) <- 1:length(rownames(gene_expression.tpm))
  colnames(gene_expression.tpm) <- c("GENEID", expr.names)
  
} else {
  # Export the gene-level expression
  gene_expression <- data.frame("GENEID"=rownames(import.gene$counts), 
                                import.gene$counts %>% as.data.frame()) %>% as.data.frame()
  rownames(gene_expression) <- 1:length(rownames(gene_expression))
  colnames(gene_expression) <- c("GENEID", expr.names)
  
  
  gene_expression.tpm <- data.frame("GENEID"=rownames(import.gene$abundance), 
                                    import.gene$abundance %>% as.data.frame()) %>% as.data.frame()
  rownames(gene_expression.tpm) <- 1:length(rownames(gene_expression.tpm))
  colnames(gene_expression.tpm) <- c("GENEID", expr.names)
}

### ----------------------------------------- Outputs -----------------------------------------
suppressWarnings(dir.create("results/tr_2_prot", recursive = T))
if (exists("snakemake")) {
  write.fasta(sequences = expr_prot_fasta,
              names = names(expr_prot_fasta), 
              file.out = unlist(snakemake@output[["expr_prot"]]))
  
  write_csv(gene_tr_prot, unlist(snakemake@output[["gene_tr_prot"]]))
  write_csv(gene_tr_prot_SP, unlist(snakemake@output[["gene_tr_prot_SP"]]))
  write_csv(expr_prot, unlist(snakemake@output[["transcript_expression"]]))
  write_csv(expr_prot_tpm, unlist(snakemake@output[["transcript_expression_tpm"]]))
  write_csv(gene_expression.tpm, unlist(snakemake@output[["gene_expression_tpm"]]))
  write_csv(gene_expression, unlist(snakemake@output[["gene_expression"]]))
} else {
  write.fasta(sequences = expr_prot_fasta,
              names = names(expr_prot_fasta), 
              file.out = "results/tr_2_prot/proteome_expressed.fasta")
  
  write_csv(gene_tr_prot, "results/tr_2_prot/gene_tr_prot.csv")
  write_csv(gene_tr_prot_SP, "results/tr_2_prot/gene_tr_prot_SP.csv")
  write_csv(expr_prot, "results/trascriptome_quant/transcript_expression.csv")
  write_csv(expr_prot_tpm, "results/trascriptome_quant/transcript_expression_tpm.csv")
  write_csv(gene_expression.tpm, "results/trascriptome_quant/gene_expression.csv")
  write_csv(gene_expression, "results/trascriptome_quant/gene_expression_tpm.csv")
}
