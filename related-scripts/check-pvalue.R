rm(list = ls())

count.tab <- read.table("/home/haoliang.xue/Documents/kamrat-test/data/test-counts.tsv", 
                        header = T, row.names = 1)

smp.sum <- colSums(count.tab)
smp.nf <- mean(smp.sum) / smp.sum

for (i in 1 : ncol(count.tab)) {
    count.tab[, i] <- count.tab[, i] * smp.nf[i]
}

cond1.cols <- c("A", "B", "C")
cond2.cols <- c("D", "E", "F", "G")

count.tab$norm.meanA <- rowMeans(count.tab[, cond1.cols])
count.tab$norm.meanB <- rowMeans(count.tab[, cond2.cols])

count.tab$p.raw <- apply(count.tab, MARGIN = 1, 
                         function(r) t.test(x = log(r[cond1.cols] + 1), 
                                            y = log(r[cond2.cols] + 1),
                                            alternative = "two.sided")$p.value)
count.tab$p.adj <- p.adjust(count.tab$p.raw, method = "BH")
