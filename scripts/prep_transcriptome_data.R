#################################################################################
# FILENAME   : prep_transcriptome_data.R
# AUTHOR     : Gulcin Ozer <ozer_gulcin@lilly.com>
# DATE       : 2025-08-11
# PROJECT    : LGM siRNA RNAseq Reports
# DESCRIPTION: Process transcriptome data and save Rdata files for each genome
#################################################################################


library(Biostrings)
library(dplyr)
library(GenomicRanges)
library(BSgenome)



## parse input arguments

args = commandArgs(trailingOnly = TRUE)

# Access arguments by index
if (length(args) > 2) {
  input.tx = args[1]
  input.name = args[2]
  out.folder = args[3]

  cat("Input taxon   :", input.tx, "\n")
  cat("Input species :", input.name, "\n")
  cat("Output folder :", out.folder, "\n")
} else {
  stop("No arguments provided. Usage: Rscript prep_transcriptome_data.R  <input_taxon> <input_species_name> <output_folder>  ")
}
#input.tx = "9606"
#input.name = "human"
#out.folder = "../data_prelim/"




# BACKGROUND INFO
# LISA transcriptome files for each genome:
# /lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/*/*_transcriptome.fasta
# Jeremy's table defining regions of transcripts
# /lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/TR_POS_TYPE.csv

# Jeremy's file that defines 3'UTRs etc
tr.pos.type = read.csv("/lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/TR_POS_TYPE.csv", sep="\t")


##### Human outlier transcript ENST00000674361.1 with length of 347561



#---------------------------------------------------------------------------
# tr.fasta
#---------------------------------------------------------------------------

#read transcriptome fasta 
tr.fasta.name = paste0("/lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/", input.tx, "/", input.tx, "_transcriptome.fasta")
tr.fasta = readDNAStringSet(tr.fasta.name)

# Split "names" into key-value pairs
parsed_list = lapply(names(tr.fasta), function(line) {
  pairs = strsplit(line, ";")[[1]]
  kv = strsplit(pairs, "=")
  keys = sapply(kv, `[`, 1)
  values = sapply(kv, `[`, 2)
  setNames(as.list(values), keys)
})

# Convert to a data frame
# details of all transcripts
tr.df = as.data.frame(bind_rows(parsed_list))
rownames(tr.df) = tr.df$transcript

write.table(tr.df, paste0(out.folder, "all_transcripts_", input.name, ".txt"), 
            sep="\t", col.names=T, row.names=F, quote=F)


# "transcript" column is the unique field for each row
# replace long annotation of fasta with the single transcript name
names(tr.fasta) = tr.df$transcript

save(tr.fasta, tr.df, file = paste0(out.folder, "transcripts_", input.name, ".RData"))



#---------------------------------------------------------------------------
# tr.fasta.UTR
#---------------------------------------------------------------------------

# Read 3'UTR definitations from Jeremy's input
dd = tr.pos.type[tr.pos.type$transcript_name %in% tr.df$transcript & tr.pos.type$transcript_pos_type=="3'UTR",]

### This takes soooo much time!!!
# Subset sequences based on coordinates
#tr.fasta.3UTR = mapply(function(name, start, end) {
#  subseq(tr.fasta[[name]], start=start, end=end)
#}, dd$transcript_name, dd$min, dd$max, SIMPLIFY =FALSE)


# getSeq from BSgenome library is much faster, and generated DNAStringSet
gr = GRanges(
  seqnames = Rle(dd$transcript_name), # transcript IDs
  ranges = IRanges(start = dd$min, end = dd$max) )

tr.fasta.UTR = getSeq(tr.fasta, gr)
#getSeq function could not maintain the seqnames, before assigning names check if widths match 
#sum(width(tr.fasta.UTR) != (dd$max - dd$min +1)) 
names(tr.fasta.UTR) = dd$transcript_name


tr.df.utr = tr.df[names(tr.fasta.UTR),]

save(tr.fasta.UTR, tr.df.utr, file = paste0(out.folder, "3pUTR_", input.name, ".RData"))


