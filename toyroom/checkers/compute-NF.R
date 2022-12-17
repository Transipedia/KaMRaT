rm(list = ls())

library(tidyr)

counts.raw <- read.table(gzfile("../data/kmer-counts.subset4toy.tsv.gz"),
                         header = TRUE)
nf.values = signif(as.numeric(1E6/colSums(counts.raw[2:ncol(counts.raw)])))
write(nf.values, "../data/nf_values.txt",
      ncolumns = length(nf.values), sep = "\t")

counts.norm <- counts.raw
for (i in 2 : ncol(counts.norm)) {
    counts.norm[, i] <- counts.norm[, i] * nf.values[i - 1]
}
counts.norm.long <- pivot_longer(counts.norm, cols = -tag,
                                 values_to = "counts", names_to = "sample")

counts.kamrat.long <- read.table("../output/recovered-kmer-counts.tsv", 
                                header = TRUE) %>%
    pivot_longer(cols = -tag,
                 values_to = "counts", names_to = "sample")

counts.cmp <- merge(counts.norm.long, counts.kamrat.long,
                    by = c("tag", "sample"))
counts.cmp$diff <- abs(counts.cmp$counts.x - counts.cmp$counts.y)
counts.cmp$relat.diff <- counts.cmp$diff / counts.cmp$counts.x
