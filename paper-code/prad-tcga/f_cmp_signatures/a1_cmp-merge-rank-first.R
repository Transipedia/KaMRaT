rm(list = ls())

library(Biostrings)
library(magrittr)

cmdArgs <- commandArgs(trailingOnly = T)
mergef.sig.path <- cmdArgs[1]
# mergef.sig.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/d_signatures/d4_rank-methods/model.train0-lrc/signatures-lasso.fa"
rankf.sig.path <- cmdArgs[2]
# rankf.sig.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/d_signatures/d4_rank-methods/model.train0-nbc/signatures-lasso.fa"
out.path <- cmdArgs[3]
# out.path <- ""

getSeqRC <- function(seq)
{
    seq.rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq.rc)
}

retrieveKMers <- function(ctg, k.len, canonical) {
    seq.len <- nchar(ctg)
    kmer.vect <- NULL
    for (start_pos in 1 : (seq.len - k.len + 1)) {
        kmer.seq <- substr(unname(ctg), start_pos, start_pos + k.len - 1)
        if (canonical) {
            kmer.seq.rc <- getSeqRC(kmer.seq)
            if (kmer.seq.rc < kmer.seq) {
                kmer.seq <- kmer.seq.rc
            }
        }
        kmer.vect <- c(kmer.vect, kmer.seq)
    }
    return(kmer.vect)
}

mergef.sig <- readDNAStringSet(mergef.sig.path)
mergef.df <- data.frame("name" = names(mergef.sig), "seq" = as.character(mergef.sig), "len" = nchar(mergef.sig))
mergef.df$mem.kmer <- lapply(mergef.df$seq, FUN = function(ctg) retrieveKMers(as.character(ctg), 31, TRUE))
# all(lapply(mergef.df$mem.kmer, FUN =length) %>% unlist() == mergef.df$len - 31 + 1)

rankf.sig <- readDNAStringSet(rankf.sig.path)
rankf.df <- data.frame("name" = names(rankf.sig), "seq" = as.character(rankf.sig), "len" = nchar(rankf.sig))
rankf.df$mem.kmer <- lapply(rankf.df$seq, FUN = function(ctg) retrieveKMers(as.character(ctg), 31, TRUE))
# all(lapply(rankf.df$mem.kmer, FUN =length) %>% unlist() == rankf.df$len - 31 + 1)

cmp.res <- c("")
for (i in 1 : nrow(mergef.df)) {
    for (j in 1 : nrow(rankf.df)) {
        if (length(intersect(unlist(mergef.df$mem.kmer[i]), unlist(rankf.df$mem.kmer[j]))) > 0) {
            cmp.res <- c(cmp.res, paste0(mergef.df$seq[i], " --- ", rankf.df$seq[j]))
        } 
    }
}
writeLines(cmp.res, paste0(out.path))

cat(paste0("nb.merge-first: ", length(mergef.sig), "\tnb.rank-first: ", length(rankf.sig), "\tnb.shared-pair: ", length(cmp.res) - 1, "\n"))
