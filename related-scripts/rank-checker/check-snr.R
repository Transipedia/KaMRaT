rm(list = ls())

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
cond.list <- as.character(unique(smp.info[, 1]))

cond1.cols <- names(count.tab)[smp.info[names(count.tab), 1] == cond.list[1]]
cond2.cols <- names(count.tab)[smp.info[names(count.tab), 1] == cond.list[2]]

if (no.norm == "no") {
    smp.sum <- colSums(count.tab)
    smp.nf <- mean(smp.sum) / smp.sum
    # write.table(smp.nf, file = paste0(out.dir, "/sample-nf.byR.txt"), col.names = F, quote = F)
} else {
    smp.nf <- rep(1, nrow(smp.info))
    names(smp.nf) <- row.names(smp.info)
}

count.tab$meanA <- apply(count.tab[, cond1.cols], MARGIN = 1, function(r) mean(r * smp.nf[cond1.cols])) %>% 
    round(digits = 6)
count.tab$meanB <- apply(count.tab[, cond2.cols], MARGIN = 1, function(r) mean(r * smp.nf[cond2.cols])) %>% 
    round(digits = 6)
count.tab$sdA <- apply(count.tab[, cond1.cols], MARGIN = 1, function(r) sd(r * smp.nf[cond1.cols])) %>% 
    round(digits = 6)
count.tab$sdB <- apply(count.tab[, cond2.cols], MARGIN = 1, function(r) sd(r * smp.nf[cond2.cols])) %>% 
    round(digits = 6)

count.tab$snr <- (count.tab$meanA - count.tab$meanB) / (count.tab$sdA + count.tab$sdB)

write.table(count.tab[order(abs(count.tab$snr), decreasing = T), ], 
	    file = paste0(out.dir, "/masked-counts.rankSNR.byR.tsv"), 
	    col.names = T, row.names = T, quote = F, sep = "\t")
