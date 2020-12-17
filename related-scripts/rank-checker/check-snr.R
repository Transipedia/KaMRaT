rm(list = ls())

library(magrittr)

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
smp.info.path <- cmdArgs[2]
out.dir <- cmdArgs[3]

# count.path <- "/home/haoliang.xue/Documents/kamrat-test/data/test-counts.tsv"
# smp.info.path <- "/home/haoliang.xue/Documents/kamrat-test/data/test-sample-conditions.tsv"

count.tab <- read.table(count.path, header = T, row.names = 1)
smp.info <- read.table(smp.info.path, header = F, row.names = 1)
cond.list <- as.character(unique(smp.info[, 1]))

cond1.cols <- names(count.tab)[smp.info[names(count.tab), 1] == cond.list[1]]
cond2.cols <- names(count.tab)[smp.info[names(count.tab), 1] == cond.list[2]]

smp.sum <- colSums(count.tab)
smp.nf <- mean(smp.sum) / smp.sum

write.table(smp.nf, file = paste0(out.dir, "/sample-nf.byR.txt"), col.names = F, quote = F)

count.tab$norm.meanA <- apply(count.tab[, cond1.cols], MARGIN = 1, 
                              function(r) mean(r * smp.nf[cond1.cols])) %>% 
    round(digits = 6)
count.tab$norm.meanB <- apply(count.tab[, cond2.cols], MARGIN = 1,
                              function(r) mean(r * smp.nf[cond2.cols])) %>% 
    round(digits = 6)
count.tab$norm.sdA <- apply(count.tab[, cond1.cols], MARGIN = 1,
                            function(r) sd(r * smp.nf[cond1.cols])) %>% 
    round(digits = 6)
count.tab$norm.sdB <- apply(count.tab[, cond2.cols], MARGIN = 1,
                            function(r) sd(r * smp.nf[cond2.cols])) %>% 
    round(digits = 6)

count.tab$snr <- (count.tab$norm.meanA - count.tab$norm.meanB) / (count.tab$norm.sdA + count.tab$norm.sdB)

write.table(count.tab[order(abs(count.tab$snr), decreasing = T), ], 
	    file = paste0(out.dir, "/masked-counts.rankSNR.byR.tsv"), 
	    col.names = T, row.names = T, quote = F, sep = "\t")
