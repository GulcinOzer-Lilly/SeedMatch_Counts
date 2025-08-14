#################################################################################
# FILENAME   : SeedMatch_counts.R
# AUTHOR     : Gulcin Ozer <ozer_gulcin@lilly.com>
# DATE       : 2025-08-12
# PROJECT    : LGM siRNA RNAseq Reports
# DESCRIPTION: Seed match counts for all possible seeds
#              Given transcriptome and seed RData variables, generate match output
#################################################################################

library(Biostrings)

## parse input arguments
args = commandArgs(trailingOnly = TRUE)

# Access arguments by index
if (length(args) > 2) {
  tr.Rdata = args[1]
  seed.Rdata = args[2]
  out = args[3]

  cat("Transcript file    :", tr.Rdata, "\n")
  cat("Seed file          :", seed.Rdata, "\n")
  cat("Output file prefix :", out, "\n")
} else {
  stop("No arguments provided. Usage: Rscript SeedMatch_counts.R <transcript Rdata> <seed Rdata> <output file prefix>  ")
}

#tr.Rdata = "../data_prelim/3pUTR_cattle.RData"
#seed.Rdata = "../data_prelim/Seeds_7mer.RData"
#out = "../data/cattle_SeedMatch_3pUTR"

#--------------------------------------------------------------------------------
# Load transcriptome
load(tr.Rdata)
load(seed.Rdata)

#rename variable for uniform processing
if( exists("dt.seed.wobble")){
   dt.seed = dt.seed.wobble
}

#rename variable for uniform processing
if( exists("tr.fasta.UTR")){
   tr.fasta = tr.fasta.UTR
   tr.df = tr.df.utr
}

seed.cols = colnames(dt.seed)[1:2]

out.by.tr = paste0(out, "_by_tr.txt")
out.by.gid = paste0(out, "_by_gid.txt")
write.table(t(c(seed.cols, 
                "taxId", "transcript", "gid", "symbol", "gene_seqname", "match.count")), 
            out.by.tr, sep="\t", col.names=F, row.names=F, quote=F)

write.table(t(c(seed.cols, "taxId", "gid", "symbol", "match.count")),
            out.by.gid, sep="\t", col.names=F, row.names=F, quote=F)


start.time <- Sys.time()


for (i in 1:dim(dt.seed)[1]){
   match.count = vcountPattern(dt.seed[i,2], tr.fasta, fixed=FALSE) 
   match.tbl = cbind(dt.seed[i,1:2], 
                     tr.df[,c("taxId", "transcript", "gid", "symbol", "gene_seqname")], 
                     match.count)
   match.tbl = match.tbl[match.tbl$match.count!=0,]
   write.table(match.tbl, out.by.tr, sep="\t", col.names=F, row.names=F, quote=F, append=T)


   match.tbl.gid = match.tbl[,c(seed.cols, "taxId", "gid", "symbol", "match.count")]
   match.tbl.gid = match.tbl.gid[order(match.tbl.gid$match.count, decreasing=T),]
   match.tbl.gid = match.tbl.gid[!duplicated(match.tbl.gid$gid),]
   write.table(match.tbl.gid, out.by.gid, sep="\t", col.names=F, row.names=F, quote=F, append=T)

}
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(paste0("Time: ", time.taken, "\n"))


# Warning messages at the end are related the rownames of tr.df
# cbind tries to use rownames, but omits due to different size 





