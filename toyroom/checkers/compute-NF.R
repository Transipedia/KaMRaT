rm(list = ls())

counts.raw <- read.table(gzfile("../data/kmer-counts.subset4toy.tsv.gz"),
                         header = TRUE, row.names = 1)
nf.values = signif(as.numeric(1E6/colSums(counts.raw)))
write(nf.values, "../data/nf_values.txt",
      ncolumns = length(nf.values), sep = "\t")
