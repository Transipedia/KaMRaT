rm (list = ls())

library(Biostrings)
library(magrittr)
library(tidyr)
library(stringr)

cmdArgs <- commandArgs(trailingOnly = T)
out.path <- cmdArgs[1]
# out.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/signature_among_train.tsv"
is_stranded <- cmdArgs[2] %>% as.logical()
# is_stranded <- TRUE
sig.path.list <- NULL
for (i in 3 : length(cmdArgs)) {
    sig.path.list <- c(sig.path.list, cmdArgs[i])
    print(paste0("fa", i - 2, ": ", cmdArgs[i]))
}
# sig.path.list <- paste0("/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train", 
#                         0 : 4, "/signatures.fa")

getSeqRC <- function(seq)
{
    seq.rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq.rc)
}

retrieveKMers <- function(contig.list, k.len, canonical) {
    tmp.df <- NULL
    for (i in 1 : length(contig.list)) {
        ctg <- contig.list[i]
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
        tmp.df <- rbind(tmp.df, 
                        data.frame("contig" = names(ctg), "kmer" = kmer.vect))
    }
    return(tmp.df)
}

kmer.contig.df <- NULL
sig.fa.list <- list()
for (i in 1 : length(sig.path.list)) {
    sig.fa <- readDNAStringSet(sig.path.list[i]) %>% as.character()
    sig.fa.list[[paste0("fa", i)]] <- sig.fa
}
kmer.contig.df <- lapply(sig.fa.list, 
                         FUN = function(s) retrieveKMers(contig.list = s, k.len = 31, canonical = !is_stranded)) %>%
    do.call(what = rbind)
kmer.contig.df$fasta.file <- str_extract(row.names(kmer.contig.df), pattern = "fa[0-9]+")
kmer.contig.df$contig_in_fasta <- paste0(kmer.contig.df$fasta.file, ":", kmer.contig.df$contig)
kmer.contig.df <- kmer.contig.df[, c("kmer", "contig_in_fasta")]
kmer.contig.df$present <- 1

sig.cmp <- aggregate(x = kmer.contig.df$contig_in_fasta, 
                     by = list(kmer.contig.df$kmer),
                     FUN = function(g) paste(sort(g), collapse = ", "))
sig.cmp <- aggregate(x = sig.cmp$Group.1, by = list(sig.cmp$x), FUN = length)
names(sig.cmp) <- c("shared.info", "shared.kmer.nb")
sig.cmp$shared.set.nb <- str_count(sig.cmp$shared.info, pattern = ",") + 1
sig.cmp <- sig.cmp[order(sig.cmp$shared.set.nb, decreasing = T), ]
write.table(sig.cmp, out.path, row.names = F, col.names = T, quote = F, sep = "\t")
