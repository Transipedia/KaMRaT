rm(list = ls())

library(Biostrings)
library(stringr)
library(magrittr)
library(compcodeR)

set.seed(91400)

cmdArgs <- commandArgs(trailingOnly = T)
nv <- cmdArgs[1] %>% as.numeric() # total feature number
es <- cmdArgs[2] %>% as.numeric() # effect size
pout <- cmdArgs[3] %>% as.numeric() # outlier percentage
out.dir <- cmdArgs[4] # output directory
# nv <- 20000
# es <- 10
# pout <- 0
# out.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-rank/compcodeR_res/"

simu.data <- generateSyntheticData(dataset = "B_625_625", n.vars = nv, effect.size = es,
                                   samples.per.cond = 150, n.diffexp = 500,
                                   fraction.upregulated = 0.5, between.group.diffdisp = FALSE,
                                   filter.threshold.total = 1, filter.threshold.mediancpm = 0,
                                   fraction.non.overdispersed = 0,
                                   random.outlier.high.prob = pout / 100, random.outlier.low.prob = 0,
                                   output.file = paste0(out.dir, "simures.rds"))
summarizeSyntheticDataSet(data.set = paste0(out.dir, "simures.rds"),
                          output.filename = paste0(out.dir, "simures.html"))

count.mat <- simu.data@count.matrix %>% as.data.frame()
count.mat$feature <- rownames(count.mat)
count.mat <- count.mat[, c("feature", paste0("sample", 1 : 300))]
v.annot <- simu.data@variable.annotations
v.annot$feature <- rownames(v.annot)
s.annot <- simu.data@sample.annotations
s.annot$sample <- rownames(s.annot)

write.table(count.mat, col.names = T, row.names = F, quote = F, sep = "\t",
            file = paste0(out.dir, "countmat.tsv"))
write.table(v.annot, col.names = T, row.names = F, quote = F, sep = "\t",
            file = paste0(out.dir, "varinfo.tsv"))
write.table(s.annot, col.names = T, row.names = F, quote = F, sep = "\t",
            file = paste0(out.dir, "smpinfo.tsv"))
