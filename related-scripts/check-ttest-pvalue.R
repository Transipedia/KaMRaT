rm(list = ls())

library(stringr)

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
out.dir <- cmdArgs[2]

# count.tab <- read.table("/home/haoliang.xue/Documents/kamrat-test/data/test-counts.tsv", 
#                         header = T, row.names = 1)

count.tab <- read.table(count.path, header = T, row.names = 1)

smp.sum <- colSums(count.tab)
smp.nf <- mean(smp.sum) / smp.sum

write.table(smp.nf, file = paste0(out.dir, "/sample-nf.byR.txt"), col.names = F, quote = F)

for (i in 1 : ncol(count.tab)) {
    count.tab[, i] <- count.tab[, i] * smp.nf[i]
}

cond1.cols <- str_detect(names(count.tab), pattern = "normal")
cond2.cols <- str_detect(names(count.tab), pattern = "tumor")

count.tab$norm.meanA <- rowMeans(count.tab[, cond1.cols])
count.tab$norm.meanB <- rowMeans(count.tab[, cond2.cols])

count.tab$p.raw <- apply(count.tab, MARGIN = 1, 
                         function(r) t.test(x = log(r[cond1.cols] + 1), 
                                            y = log(r[cond2.cols] + 1),
                                            alternative = "two.sided")$p.value)
count.tab$p.adj <- p.adjust(count.tab$p.raw, method = "BH")

write.table(count.tab[order(count.tab$p.adj, decreasing = F), ], 
	    file = paste0(out.dir, "/masked-counts.rank.byR.tsv"), 
	    col.names = T, row.names = T, quote = F, sep = "\t")
