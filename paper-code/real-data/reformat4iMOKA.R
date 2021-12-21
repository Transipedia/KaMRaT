rm (list = ls())

library(stringr)

cmdArgs <- commandArgs(trailingOnly = T)
ranked.tab.path <- cmdArgs[1]
smp.condi.path <- cmdArgs[2]
out.path <- cmdArgs[3]

# ranked.tab.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/kamrat_res/fold0/mergingfirst/merge-rank-res.pearson20.ttest.padj.tsv"
# smp.condi.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/matrices/sampleshuf.train0.tsv"
# out.path <- "/home/haoliang.xue/Downloads/tmp.tsv"

smp.condi <- read.table(smp.condi.path, header = F, row.names = 1)
colnames(smp.condi) <- "group"

ranked.tab <- read.table(ranked.tab.path, header = T, row.names = 1)
ranked.tab <- ranked.tab[, rownames(smp.condi)]

smp.condi <- as.data.frame(t(smp.condi))
smp.condi <- smp.condi[, colnames(ranked.tab)]

cat("\t", file = out.path, append = F)
write.table(smp.condi, out.path, append = T, row.names = T, col.names = T, quote = F, sep = "\t")
write.table(ranked.tab, out.path, append = T, row.names = T, col.names = F, quote = F, sep = "\t")
