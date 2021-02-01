rm(list = ls())

library(glmnet)
library(Matrix)
library(magrittr)

cmdArgs <- commandArgs(trailingOnly = T)
trained.mdl.path <- cmdArgs[1]
# trained.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0/trained.model.rds"
sig.count.path <- cmdArgs[2]
# sig.count.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/count-in-tests/contig-counts.test0.tsv"
smp.info.path <- cmdArgs[3]
# smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.test0.tsv"

trained.mdl <- readRDS(trained.mdl.path)

sig.count <- read.table(sig.count.path, header = T, row.names = 1) %>%
    t() %>%
    as.data.frame()
sig.count <- sig.count[, rownames(coefficients(trained.mdl))[-1]]

smp.info <- read.table(smp.info.path, header = T, row.names = 1)
smp.info <- smp.info[row.names(sig.count), ]

auc.pred <- assess.glmnet(trained.mdl, 
                          newx = data.matrix(sig.count), newy = smp.info$condition,
                          family = "binomial", s = "lambda.1se")$auc
print(paste0("AUC of ", trained.mdl.path, ": ", round(auc.pred, 3)))
