rm(list = ls())

library(stringr)

kmer.tab.sel.path <- "/home/haoliang.xue/development/kamrat-test/data/luad-diff-counts-up-sel.tsv"
merged.contig.path <- "/home/haoliang.xue/development/kamrat-test/merged.counts.tsv"

kmer.tab.sel <- read.table(kmer.tab.sel.path, header = T)
merged.contig <- read.table(merged.contig.path, header = T)

contig.counts <- as.numeric(merged.contig[, str_detect(names(merged.contig), "ERR")])
kmer.mean.counts <- as.numeric(apply(kmer.tab.sel[, str_detect(names(kmer.tab.sel), "ERR")], MARGIN = 2, mean))

d <- abs(contig.counts - kmer.mean.counts)
max(d)
rd <- d
for (i in 1 : length(rd)) {
    if (kmer.mean.counts[i] != 0) {
        rd[i] <- rd[i] / kmer.mean.counts[i]
    }
}
max(rd)
which.max(rd)
contig.counts[79]
kmer.mean.counts[79]
