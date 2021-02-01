rm(list = ls())

library(stringr)
library(magrittr)

source("/store/USERS/haoliang.xue/development/KaMRaT/related-scripts/rank-checker/calc-scores.func.R")

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
# count.path <- "/home/haoliang.xue/development/kamrat-test/data/test-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/development/kamrat-test/data/test-sample-conditions.tsv"
out.dir <- cmdArgs[3]
# out.dir <- ""
no.norm <- cmdArgs[4] %>% as.logical()
# no.norm <- "no"

smp.info <- read.table(smp.info.path, header = F, row.names = 1)
count.tab <- read.table(count.path, header = T, row.names = 1)[, rownames(smp.info)]

if (!no.norm) {
    smp.sum <- colSums(count.tab)
    smp.nf <- mean(smp.sum) / smp.sum
} else {
    smp.nf <- rep(1, nrow(smp.info))
    names(smp.nf) <- row.names(smp.info)
}

feature.score <- data.frame(row.names = rownames(count.tab), 
                            "relat_sd" = calcRelatSD(count.tab, smp.nf),
                            "ttest.padj" = calcTtestPvalue(count.tab, smp.info, smp.nf),
                            "snr" = calcSNR(count.tab, smp.info, smp.nf),
                            "lrc.f1" = calcLRCF1(count.tab, smp.info, smp.nf),
                            "nbc.f1" = calcNBCF1(count.tab, smp.info, smp.nf),
                            "svm.hingeloss" = calcSVMHingeLoss(count.tab, smp.info, smp.nf))

write.table(feature.score, file = paste0(out.dir, "/feature-scores.byR.tsv"), 
            col.names = T, row.names = T, quote = F, sep = "\t")
