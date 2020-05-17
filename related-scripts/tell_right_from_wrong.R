rm(list = ls())

library(Biostrings)

cmdArgs <- commandArgs(trailingOnly = T)

align.tsv <- cmdArgs[1]
seq.fa <- cmdArgs[2]
out.fa <- cmdArgs[3]

align.info <- read.table(align.tsv, header = T)
contig.ori <- readDNAStringSet(seq.fa)

contig.right.name <- align.info$qseqid[align.info$qstart == 1 & 
                                       align.info$qend == align.info$qlen &
                                       align.info$nident == align.info$qlen &
                                       align.info$pident == 100.000]

contig.right <- contig.ori[names(contig.ori) %in% contig.right.name]
names(contig.right) <- paste0(names(contig.right), ":right")
contig.wrong <- contig.ori[!names(contig.ori) %in% contig.right.name]
names(contig.wrong) <- paste0(names(contig.wrong), ":wrong")

contig.final <- c(contig.right, contig.wrong)
writeXStringSet(contig.final, out.fa, width = 20001)