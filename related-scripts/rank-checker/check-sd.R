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

cols <- c("A", "B", "C", "D", "E", "F", "G")

count.tab$sd <- apply(count.tab[, cols], MARGIN = 1, sd)
count.tab$rsd <- apply(count.tab[, cols], MARGIN = 1, function(r) sd(r)/mean(r))

write.table(count.tab[order(abs(count.tab$rsd), decreasing = T), ], 
            file = paste0(out.dir, "/masked-counts.rankSD.byR.tsv"), 
            col.names = T, row.names = T, quote = F, sep = "\t")
