rm(list = ls())

library(stringr)

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
out.dir <- cmdArgs[2]

count.path <- "/home/haoliang.xue/Documents/kamrat-test/data/test-counts.tsv"

count.tab <- read.table(count.path, header = T, row.names = 1)

smp.sum <- colSums(count.tab)
smp.nf <- mean(smp.sum) / smp.sum

for (i in 1 : ncol(count.tab)) {
    count.tab[, i] <- count.tab[, i] * smp.nf[i]
}

# cond1.cols <- str_detect(names(count.tab), pattern = "normal")
# cond2.cols <- str_detect(names(count.tab), pattern = "tumor")

cond1.cols <- c("A", "B", "C")
cond2.cols <- c("D", "E", "F", "G")

count.tab$norm.meanA <- rowMeans(count.tab[, cond1.cols])
count.tab$norm.meanB <- rowMeans(count.tab[, cond2.cols])
count.tab$norm.sdA <- apply(count.tab[, cond1.cols], MARGIN = 1, sd)
count.tab$norm.sdB <- apply(count.tab[, cond2.cols], MARGIN = 1, sd)

count.tab$snr <- (count.tab$norm.meanA - count.tab$norm.meanB) / (count.tab$norm.sdA + count.tab$norm.sdB)

write.table(count.tab[order(abs(count.tab$snr), decreasing = T), ], 
	    file = paste0(out.dir, "/masked-counts.rank.byR.tsv"), 
	    col.names = T, row.names = T, quote = F, sep = "\t")
