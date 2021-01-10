rm(list = ls())

library(stringr)
library(magrittr)

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
smp.info.path <- cmdArgs[2]
out.dir <- cmdArgs[3]
no.norm <- cmdArgs[4]

# count.path <- "/home/haoliang.xue/development/kamrat-test/data/test-counts.tsv"
# smp.info.path <- "/home/haoliang.xue/development/kamrat-test/data/test-sample-conditions.tsv"
# no.norm <- TRUE

count.tab <- read.table(count.path, header = T, row.names = 1)
smp.info <- read.table(smp.info.path, header = F, row.names = 1)

if (no.norm == "no") {
    smp.sum <- colSums(count.tab)
    smp.nf <- mean(smp.sum) / smp.sum
    # write.table(smp.nf, file = paste0(out.dir, "/sample-nf.byR.txt"), col.names = F, quote = F)
} else {
    smp.nf <- rep(1, nrow(smp.info))
    names(smp.nf) <- row.names(smp.info)
}

count.tab$mean <- apply(count.tab[, names(smp.nf)], MARGIN = 1, function(r) mean(r * smp.nf))
count.tab$sd <- apply(count.tab[, names(smp.nf)], MARGIN = 1, function(r) sd(r * smp.nf))
count.tab$rsd <- count.tab$sd / count.tab$mean

write.table(count.tab[order(abs(count.tab$rsd), decreasing = T), ],
	    file = paste0(out.dir, "/masked-counts.rankRSD.byR.tsv"), 
            col.names = T, row.names = T, quote = F, sep = "\t")
