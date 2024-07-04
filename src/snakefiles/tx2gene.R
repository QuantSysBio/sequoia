### ---------------------------------------------- Transcript-gene mapping  ----------------------------------------------
# description:  Create a gene to transcript mapping including stringtie de-novo assembly transcripts and loci.
#               Reference annotation is used where possible. 
#               Novel isoforms of known genes have new transcript IDs and reference gene IDs.
#               New gene IDs only used for novel loci.
#
# input:        1. Reference transcriptome (.fasta and .gtf)
#               2. Gffcompare transcript tracking table
# output:       
#               - Transcript to gene mapping table.
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(dplyr)
library(stringr)
library(rtracklayer)

# Load the trasfrag tracking from gffcompare
gffcomp_tracking <- read.delim(snakemake@input[["gffcomp_tracking"]], header = F)
colnames(gffcomp_tracking) <- c("transcript","locus", "ref_gene", "class_code", "ref_tr")

# Load the transcriptome fasta
tr <- seqinr::read.fasta(snakemake@input[["transcripts_fa"]])

tr_tcons <- grep(pattern = "TCONS_", x = names(tr), fixed = T)
tr_ref <- tr[-tr_tcons]
tr_tcons <- tr[tr_tcons]

# Check if all the novel transcripts are present in the tracking file:
setdiff(names(tr_tcons), gffcomp_tracking$transcript)
track_tcons <- gffcomp_tracking[gffcomp_tracking$transcript %in% names(tr_tcons),]


# ------------------------- Transript/Gene matching for tximport  -------

# 1. If the transcript is in the reference, use the reference gene/tr match
ref_gtf <- rtracklayer::import(snakemake@input[["gencode_annotation_gtf"]]) %>% as_tibble()
tx2gene <- data.frame("TXNAME"=ref_gtf$transcript_id,
                      "GENEID"=ref_gtf$gene_id) %>% unique() %>% na.omit()

# 2. If the transcript is novel, but has an annotated reference transcript - use the reference parent gene
tr_names <- track_tcons$ref_gene %>% str_split(pattern = "\\|", n=Inf,simplify = T) %>% as.data.frame()
tr_names <- cbind(track_tcons$transcript, tr_names)

tr_match <- tr_names[!tr_names$V1 == "-",]
tr_novel <- tr_names[tr_names$V1 == "-",]

colnames(tr_match)[3] <- "TXNAME"
tr_match <- left_join(tr_match, tx2gene)
tr_match <- tr_match[,c(1,4)]
colnames(tr_match)[1] <- "TXNAME"

# 3. Use MSTRG gene IDs from stringtie.merge for the remaining transcripts
tr_tcons_novel <- track_tcons[track_tcons$transcript %in% tr_novel$`track_tcons$transcript`,]
tr_tcons_novel_reftr <- tr_tcons_novel$ref_tr %>% str_split(pattern = "\\|", n=Inf,simplify = T) %>% as.data.frame()
tr_tcons_novel_reftr$V1 <- tr_tcons_novel_reftr$V1 %>% gsub(pattern = "q1:", replacement = "", fixed = T)

tr_novel <- data.frame("TXNAME"=tr_novel$`track_tcons$transcript`,
                       "GENEID"=tr_tcons_novel_reftr$V1) %>% unique() %>% na.omit

# Join three reference tables together
tx2gene <- rbind(tx2gene, tr_match, tr_novel)

# ------------------------- 
# Export results
write.csv(tx2gene, file = unlist(snakemake@output[["tx2gene"]]))
