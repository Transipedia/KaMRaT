rm(list = ls())

library(glmnet)
library(Matrix)
library(caret)

cmdArgs <- commandArgs(trailingOnly = T)
trained.mdl.path <- cmdArgs[1]
trained.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0/trained.model.rds"
sig.count.path <- cmdArgs[2]
sig.count.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/count-in-tests/contig-counts.test0.tsv"
smp.info.path <- cmdArgs[3]
smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.test0.tsv"

trained.mdl <- readRDS(trained.mdl.path)
sig.count <- read.table(sig.count.path, header = T, row.names = 1)
smp.info <- read.table(smp.info.path, header = T)

print(paste("signature matrix's order same as that in model:", all(rownames(sig.count) == rownames(trained.mdl$beta))))

mdl.eval <- predict(mdl.postsel,
                    newx = sig.count,
                    family = "bionomial",
                    s = mdl.postsel$lambda.1se)
mdl.eval$auc
