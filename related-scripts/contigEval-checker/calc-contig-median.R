rm(list = ls())

library(Biostrings)
library(Matrix)
library(magrittr)

cmdArgs <- commandArgs(trailingOnly = T)
contig.path <- cmdArgs[1]
# contig.path <- "/home/haoliang.xue/development/kamrat-test/data/test-contig.fasta"
kmer.mat.path <- cmdArgs[2]
# kmer.mat.path <- "/home/haoliang.xue/development/kamrat-test/data/test-counts.tsv"
output.path <- cmdArgs[3]
stranded <- cmdArgs[4] %>% as.logical()
# stranded <- TRUE
sample.info.path <- cmdArgs[5]
# sample.info.path <- "/home/haoliang.xue/development/kamrat-test/data/test-sample-conditions.tsv"

countContigMedian <- function(contig.seq, kmer.mat, k.len) {
    mem.kmer.mat <- NULL
    for (start_pos in 1 : (nchar(contig.seq) - k.len + 1)) {
        kmer <- substr(contig.seq, start_pos, start_pos + k.len - 1)
        if (kmer %in% row.names(kmer.mat)) {
            mem.kmer.mat <- rbind(mem.kmer.mat, kmer.mat[kmer, ])
        }
    }
    if (is.null(mem.kmer.mat)) {
        contig.median <- matrix(rep(0, ncol(kmer.mat)), nrow = 1) %>% 
            as.data.frame()
        names(contig.median) <- names(kmer.mat)
    } else {
        contig.median <- apply(mem.kmer.mat, MARGIN = 2, median) %>% 
            t() %>% 
            as.data.frame()
    }
    contig.median$contig <- unname(contig.seq)
    return(contig.median[, c("contig", names(mem.kmer.mat))])
}

contigs <- readDNAStringSet(contig.path) %>% as.character()
kmer.mat <- read.table(kmer.mat.path, header = T, row.names = 1)
smp.info <- read.table(sample.info.path, header = F, row.names = 1)
kmer.mat <- kmer.mat[, colnames(kmer.mat) %in% rownames(smp.info)]
nf.vect <- colSums(kmer.mat)
nf.vect <- mean(nf.vect) / nf.vect
for (j in 1 : ncol(kmer.mat)) {
    kmer.mat[, j] <- kmer.mat[, j] * nf.vect[j]
}
ctg.count <- lapply(contigs, FUN = function(c) countContigMedian(contig.seq = c, kmer.mat = kmer.mat, k.len = 31)) %>%
    do.call(what = rbind)
write.table(ctg.count, file = output.path, row.names = F, col.names = T, quote = F, sep = "\t")