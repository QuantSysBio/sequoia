### SPI-art project ###
# 
# description:  In silico translate transcriptome UTRs, filter by expression, translate in 3 frames
# input:        - human Gencode transcriptome annotation .gtf
#               - human Gencode genome package 
# output:       A folder with: .fasta of 3-frame translation.
#               
# author:       Yehor Horokhovskyi

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(tidyr)
library(ORFik)  
library(plyranges)
library(stringr)
library(vroom)

###-------------------Creating the result folder--------------------------
create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created:", folder_path, "\n")
  } else {
    cat("Folder already exists:", folder_path, "\n")
  }
}
folder_path <- 'results/Exhaustive/cryptic/all/'
create_folder_if_not_exists(folder_path)

### ------------------ Start user variables ------------------------------

# Multi-threading
Ncpu = 7

# For non-CDS ORF discovery
min_ORF_length = 8
alt_init_codons = snakemake@params[['init_codons']]
# attributes(GENETIC_CODE)$alt_init_codons <- c("TTG", "CTG") # original
attributes(GENETIC_CODE)$alt_init_codons <- strsplit(snakemake@params[['init_codons']], '|')[[1]] # cryptic

# GTF annotation
input <- list()
input$gtf_all <- read_gff(file = paste0(snakemake@input[['reference_gff3']]),
                          genome_info = "hg38")
input$gtf_all

# Genome object
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
Hsapiens

# update the names to match
upd_names <- names(Hsapiens) %>% 
  str_replace_all(pattern = "v1", replacement = "\\.1") %>% 
  str_replace_all(pattern = "v2", replacement = "\\.2") %>% 
  str_extract(pattern = str_c(unique(as.character(seqnames(input$gtf_all)))[26:length(unique(as.character(seqnames(input$gtf_all))))], collapse = "|"))
seqnames(Hsapiens)[which(!is.na(upd_names))] <- upd_names[which(!is.na(upd_names))]

# txdb <- makeTxDbFromGRanges(input$gtf_all)
txdb <- makeTxDbFromGFF(file = paste0(snakemake@input[['reference_gff3']]),
                        dataSource="GENCODE v45",
                        organism="Homo sapiens")
txdb

### ------------------ End user variables ------------------------------
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Extract ORFs from genomic features
seqlevels_txdb <- seqlevels(txdb)
seqlevels_txdb <- seqlevels_txdb[!seqlevels_txdb=="chrM"]
seqlevels_txdb <- seqlevels_txdb[str_starts(seqlevels_txdb, "chr")]
seqlevels_txdb

genomic_strata <- c("fiveUTR", "threeUTR", "lncRNA", "intronic", "intergenic")

for (stratum in genomic_strata[1]) {
  print(stratum)
  for (chromosome in seqlevels_txdb[1]) {
    # Process one chromosome at the time
    print(chromosome)
    print(Sys.time())
    seqlevels(txdb) <- chromosome
    
    # Extract CDS coordinates
    CDS <- cdsBy(txdb, by = "tx", use.names=TRUE)
    CDS_seq <- extractTranscriptSeqs(Hsapiens, CDS)
    
    ### ------------------ Extract and map ORFs ------------------------------
    switch(stratum,
           fiveUTR={
             # Extract 5'-UTR coordinates
             stratum_coords <- fiveUTRsByTranscript(txdb) %>%
               as_granges() %>%
               select(-c("exon_rank", "group", "group_name", "exon_id")) %>%
               unique() %>%
               group_by(exon_name) %>%
               reduce_ranges_directed()
             stratum_coords <- split(stratum_coords, stratum_coords$exon_name) 
             stratum_seq <- extractTranscriptSeqs(Hsapiens, stratum_coords)
           },
           threeUTR={
             # Extract 3'-UTR coordinates
             stratum_coords <- threeUTRsByTranscript(txdb) %>%
               as_granges() %>%
               select(-c("exon_rank", "group", "group_name", "exon_id")) %>%
               unique() %>%
               group_by(exon_name) %>%
               reduce_ranges_directed()
             stratum_coords <- split(stratum_coords, stratum_coords$exon_name) 
             stratum_seq <- extractTranscriptSeqs(Hsapiens, stratum_coords)
           },
           lncRNA={
             # Extract lncRNA coordinates
             stratum_coords <- input$gtf_all %>% 
               filter(transcript_type == "lncRNA" & type == "transcript" & seqnames == chromosome) %>% 
               as_granges() %>%
               select(transcript_id) %>%
               unique()
             stratum_coords <- split(stratum_coords, stratum_coords$transcript_id) 
             stratum_seq <- extractTranscriptSeqs(Hsapiens, stratum_coords)
           },
           intronic={
             # Extract intron coordinates
             stratum_coords <- intronsByTranscript(txdb, use.names=TRUE) %>%
               as_granges() %>%
               select(group_name) %>%
               unique() %>%
               group_by(group_name) %>%
               reduce_ranges_directed() %>% 
               as_tibble() %>%
               group_by(group_name) %>%
               mutate(group_name = paste0(group_name, "_", row_number())) %>%
               as_granges()
             stratum_coords <- split(stratum_coords, stratum_coords$group_name) 
             stratum_seq <- extractTranscriptSeqs(Hsapiens, stratum_coords)
           },
           intergenic={
             gene_regions = exonsBy(txdb, "gene")
             stratum_coords = gaps(unlist(range(gene_regions))) 
             # %>%
             #   group_by(seqnames, end, strand) %>%
             #   reduce_ranges_directed() 
             
             stratum_seq <- getSeq(Hsapiens, stratum_coords)
             names(stratum_seq) <- paste(as.character(stratum_coords@seqnames),
                                         as.character(stratum_coords@ranges),
                                         as.character(stratum_coords@strand), sep = "_")
             
             names(stratum_coords) <- names(stratum_seq)
             stratum_coords <- as(stratum_coords, "GRangesList")
           },
           stop("No method defined for a given stratum")
    )
    
    if (length(stratum_coords) > 0) {
      ### ------------------ Find ORFs ------------------------------
      stratum_ORFs <- findMapORFs(grl = stratum_coords, 
                                  seqs = stratum_seq, 
                                  startCodon = alt_init_codons, 
                                  longestORF = T, 
                                  minimumLength = min_ORF_length, 
                                  groupByTx = FALSE) %>%
        as_granges() %>%
        mutate(group_name = str_split_fixed(group_name, "_", 2)[,1]) %>%
        mutate(group = as.integer(as.factor(group_name))) 
      
      # 1. Which ORFs overlap CDS
      stratum_ORFs_CDS <- find_overlaps_directed(stratum_ORFs, as_granges(CDS), suffix = c(".x", ".y"))
      
      ### 2. Whether an ORF is in-frame with CDS
      CDS_tbl <- CDS %>%
        as_tibble() %>%
        rename(group_name.y = group_name,
               start.y = start,
               end.y = end) %>%
        select(group_name.y, start.y, end.y, cds_id, strand, cds_id, cds_name, exon_rank) 
      
      stratum_ORFs_CDS <- stratum_ORFs_CDS %>% 
        as_tibble() %>%
        left_join(CDS_tbl) %>%
        mutate(start_relative = (start - start.y))  %>%
        mutate(end_relative = (end - end.y)) %>%
        mutate(CDS_in_frame = is.wholenumber(end_relative/3)) %>%
        as_granges()
      
      stratum_ORFs_CDS_inframe <- stratum_ORFs_CDS %>%
        filter(CDS_in_frame == T)  %>% 
        filter(!(end_relative == 0 & start_relative == 0)) %>%
        select(-exon_rank)
      
      ### 3. Annotate overlaps with CDS
      stratum_ORF_overlap_coords <- stratum_ORFs[stratum_ORFs$names %in% stratum_ORFs_CDS$names] %>%
        join_overlap_intersect_directed(as_granges(CDS)) %>%
        as_tibble() %>%
        filter(names %in% stratum_ORFs_CDS_inframe$names) %>%
        select(seqnames, start, end, strand, names, group_name.y) %>%
        rename(start.overlap = start, 
               end.overlap = end)
      
      stratum_ORFs_CDS_inframe <- stratum_ORFs_CDS_inframe  %>%
        as_tibble() %>%
        left_join(stratum_ORF_overlap_coords) %>%
        select(-c("cds_id", "cds_name", "group.y")) %>%
        as_granges()
      
      ### ------------------ Prepare outputs ------------------------------
      ### 4. Extract sequences
      annotation_tables <- list()
      fasta_files <- list()
      
      # No CDS overlap
      stratum_ORFs_noCDS <- stratum_ORFs[!stratum_ORFs$names %in% stratum_ORFs_CDS$names] %>%
        mutate(stratum = stratum) %>%
        mutate(overlaps_CDS = "no_CDS_overlaps") %>%
        as_tibble() %>%
        select(-c("group", "group_name")) %>% 
        filter(width >= 3*min_ORF_length)  %>%
        # mutate(expressed = ifelse(str_split_fixed(names, "_", 2)[,1] %in% expr_tr_exon, TRUE, FALSE)) %>%
        mutate(ORF_name = paste(names,seqnames,start,end,strand,stratum,overlaps_CDS, sep = "|"))
      annotation_tables$ORFs_noCDS <- stratum_ORFs_noCDS
      
      fasta_files$ORFs_noCDS <- stratum_ORFs_noCDS %>%
        # select(seqnames, start, end, width, strand, ORF_name, expressed) %>%
        select(seqnames, start, end, width, strand, ORF_name) %>%
        unique() %>%
        as_granges()
      
      fasta_files$ORFs_noCDS <- fasta_files$ORFs_noCDS %>%
        split(fasta_files$ORFs_noCDS$ORF_name) %>%
        extractTranscriptSeqs(x = Hsapiens)
      
      # CDS off-frame
      annotation_tables$ORFs_CDS_frameshift <- stratum_ORFs_CDS %>%
        filter(CDS_in_frame == F) %>%
        select(-c("group.x", "group_name.x", "exon_rank", "cds_id", "cds_name", "CDS_in_frame", "group.y")) %>%
        as_tibble() %>%
        rename(start.CDS = start.y,
               end.CDS = end.y,
               transcript.CDS = group_name.y) %>%
        mutate(stratum = stratum) %>%
        mutate(overlaps_CDS = "out-of-frame") %>% 
        filter(width >= 3*min_ORF_length) %>%
        # mutate(expressed = ifelse(str_split_fixed(names, "_", 2)[,1] %in% expr_tr_exon, TRUE, FALSE)) %>%
        mutate(ORF_name = paste(names,seqnames,start,end,strand,stratum,overlaps_CDS, sep = "|"))
      
      fasta_files$ORFs_CDS_frameshift <- annotation_tables$ORFs_CDS_frameshift %>%
        # select(seqnames, start, end, width, strand, ORF_name, expressed) %>%
        select(seqnames, start, end, width, strand, ORF_name) %>%
        unique() %>%
        as_granges()
      
      fasta_files$ORFs_CDS_frameshift <- fasta_files$ORFs_CDS_frameshift %>%
        split(fasta_files$ORFs_CDS_frameshift$ORF_name) %>%
        extractTranscriptSeqs(x = Hsapiens)
      
      # CDS in-frame
      annotation_tables$ORFs_CDS_inframe <- stratum_ORFs_CDS_inframe %>%
        as_tibble() %>%
        select(-c("group.x", "group_name.x", "CDS_in_frame")) %>%
        rename(start.CDS = start.y,
               end.CDS = end.y,
               transcript.CDS = group_name.y) %>%
        mutate(stratum = stratum) %>%
        mutate(overlaps_CDS = "in-frame") %>% 
        filter(width >= 3*min_ORF_length)  %>%
        # mutate(expressed = ifelse(str_split_fixed(names, "_", 2)[,1] %in% expr_tr_exon, TRUE, FALSE)) %>%
        mutate(ORF_name = paste(names,seqnames,start,end,strand,stratum,overlaps_CDS, sep = "|"))
      
      fasta_files$ORFs_CDS_inframe <- annotation_tables$ORFs_CDS_inframe  %>%
        # select(seqnames, start, end, width, strand, ORF_name, expressed) %>%
        select(seqnames, start, end, width, strand, ORF_name) %>%
        unique() %>%
        as_granges()
      
      fasta_files$ORFs_CDS_inframe <- fasta_files$ORFs_CDS_inframe %>%
        split(fasta_files$ORFs_CDS_inframe$ORF_name) %>%
        extractTranscriptSeqs(x = Hsapiens)
      
      ### Export FASTA
      # All - AA
      fasta_files$ORFs_noCDS %>%
        translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_noCDS.fa"), append = T)
      
      fasta_files$ORFs_CDS_frameshift %>%
        translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_frameshift.fa"), append = T)
      
      fasta_files$ORFs_CDS_inframe %>%
        translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_inframe.fa"), append = T)
      
      fasta_files$ORFs_noCDS %>%
        # translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_noCDS_NT.fa"), append = T)
      
      # All - NT
      fasta_files$ORFs_CDS_frameshift %>%
        # translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_frameshift_NT.fa"), append = T)
      
      fasta_files$ORFs_CDS_inframe %>%
        # translate(no.init.codon = T, if.fuzzy.codon = "solve") %>%
        unique() %>%
        writeXStringSet(filepath = paste0(snakemake@output[['all_result']], stratum, "_ORFs_inframe_NT.fa"), append = T)
      
      ### Export tables
      # All
      annotation_tables$ORFs_noCDS %>%
        select(-ORF_name) %>%
        unique() %>%
        vroom_write(file = paste0(snakemake@output[['all_result']], stratum, "_ORFs_noCDS.csv"), delim = ",", append = T, num_threads = Ncpu)
      
      annotation_tables$ORFs_CDS_frameshift %>%
        select(-ORF_name) %>%
        unique() %>%
        vroom_write(file = paste0(snakemake@output[['all_result']], stratum, "_ORFs_frameshift.csv"), delim = ",", append = T, num_threads = Ncpu)
      
      annotation_tables$ORFs_CDS_inframe %>%
        select(-ORF_name) %>%
        unique() %>%
        vroom_write(file = paste0(snakemake@output[['all_result']], stratum, "_ORFs_inframe.csv"), delim = ",", append = T, num_threads = Ncpu)
    }
  }
}

### ------------------ Add colnames to output tables ------------------------------
dir_res <- c(snakemake@output[['all_result']]#, "results/expressed/"
             )
for (dir in dir_res) {
  
  ### Remove empty files
  list.files(dir,  full.names = T)[file.size(list.files(dir,  full.names = T)) == 0] %>%
    file.remove()
  
  output_tables <- list.files(dir, pattern = ".csv")
  for (i in seq_along(output_tables)) {
    df <- vroom(paste0(dir, output_tables[[i]]), col_names = F, num_threads = Ncpu, show_col_types = FALSE)
    
    if (nrow(df) > 0) {
      switch(as.character(ncol(df)),
             "8"={
               colnames(df) <- annotation_tables$ORFs_noCDS %>%
                 select(-ORF_name) %>%
                 colnames()
             },
             "13"={
               colnames(df) <- annotation_tables$ORFs_CDS_frameshift %>%
                 select(-ORF_name) %>%
                 colnames()
             },
             "15"={
               colnames(df) <- annotation_tables$ORFs_CDS_inframe %>%
                 select(-ORF_name) %>%
                 colnames()
             },
             stop("No method defined for a given column number")
      )
      vroom_write(df, paste0(dir, output_tables[[i]]), delim = ",", col_names = TRUE, num_threads = Ncpu)
      
    } else {
      file.remove(paste0(dir, output_tables[[i]]))
    }
    rm(df)
  }
}

file.create(snakemake@output[['Workflow_output']])
