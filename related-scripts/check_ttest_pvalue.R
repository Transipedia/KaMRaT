rm(list = ls())

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
out.path <- cmdArgs[2]

# count.path <- "/home/haoliang.xue/Documents/kamrat-test/data/test-counts.tsv"
out.path <- "/home/haoliang.xue/Documents/kamrat-test/kamrat-rank/tmp.tsv"

count.tab <- read.table(count.path, header = T, row.names = 1)
count.tab <- count.tab + 1
smp.sum <- colSums(count.tab)
smp.sum.mean <- mean(smp.sum)
for (i in 1 : ncol(count.tab)) {
    count.tab[, i] <- count.tab[, i] / smp.sum[i] * smp.sum.mean
}
count.tab$norm.meanA <- rowMeans(count.tab[, c("A", "B", "C")])
count.tab$norm.meanB <- rowMeans(count.tab[, c("D", "E", "F", "G")])

count.tab$score.byR <- apply(count.tab, MARGIN = 1,
                             function(r) t.test(x = log(r[c("A", "B", "C")]), 
                                                y = log(r[c("D", "E", "F", "G")]),
                                                alternative = "two.sided")$p.value)
count.tab$score.byR <- p.adjust(count.tab$score.byR, method = "BH")

write.table(count.tab[order(count.tab$score, decreasing = F), ], 
            file = out.path, col.names = T, row.names = T, quote = F, sep = "\t")
