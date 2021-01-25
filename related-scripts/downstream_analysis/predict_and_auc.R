rm(list = ls())

library(glmnet)
library(Matrix)
library(caret)

cmdArgs <- commandArgs(trailingOnly = T)
trained.mdl.path <- cmdArgs[1]
# trained.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/c_kamrat_contigs/model.train0/trained.model.rds"
sig.count.path <- cmdArgs[2]
# sig.count.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/c_kamrat_contigs/count-in-tests/contig-counts.test0.tsv"
smp.info.path <- cmdArgs[3]
# smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.test0.tsv"

trained.mdl <- readRDS(trained.mdl.path)
sig.count <- read.table(sig.count.path, header = T, row.names = 1) %>%
    t() %>% 
    data.matrix()
smp.info <- read.table(smp.info.path, header = T)
smp.info <- smp.info[, c("sample", "condition")]

print(paste("signature matrix row names same order as that in model:", all(names(sig.count) == rownames(trained.mdl$beta))))
print(paste("signature matrix row names same order as sample info:", all(names(sig.count) == smp.info$condition)))

mdl.eval <- assess.glmnet(trained.mdl[["trained.mdl"]], newx = sig.count, newy = smp.info$condition, 
                          family = "binomial", s = trained.mdl[["lambda.1se"]])
print(paste0("AUC of ", trained.mdl.path, ": ", as.numeric(mdl.eval$auc)))
