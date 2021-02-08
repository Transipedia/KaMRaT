rm(list = ls())

library(stringr)
library(Matrix)
library(magrittr)
library(glmnet)
library(randomForest)
library(ROCR)

cmdArgs <- commandArgs(trailingOnly = T)
trained.mdl.path <- cmdArgs[1]
# trained.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0/fitted-model-randomforest.rds"
sig.count.path <- cmdArgs[2]
# sig.count.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/count-in-tests/contig-counts.test0-randomforest.tsv"
smp.info.path <- cmdArgs[3]
# smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.test0.tsv"

evalModel <- function(mdl.fit, test.x, test.y, mdl.name) {
    if (mdl.name == "lasso" || mdl.name == "elasticnet" || mdl.name == "ridge") {
        pred <- predict(mdl.fit, newx = test.x, type = "response")
        pred4auc <- prediction(predictions = pred, labels = test.y)
        auc.pred <- performance(pred4auc, measure = "auc")@y.values[[1]]
    } else if (mdl.name == "randomforest") {
        pred <- predict(mdl.fit, test.x, type = "prob")
        pred4auc <- prediction(predictions = pred[, 2], 
                               labels = ifelse(test.y == colnames(pred)[2], yes = 1, no = 0))
        auc.pred <- performance(pred4auc, measure = "auc")@y.values[[1]]
    } else {
        stop(paste0("unknown model name to fit: ", mdl.name, " (acceptable options: lasso, elasticnet, ridge, randomforest)"))
    }
    return(auc.pred)
}

model.name <- str_extract(trained.mdl.path, pattern = "lasso|elasticnet|ridge|randomforest")

trained.mdl <- readRDS(trained.mdl.path)

smp.info <- read.table(smp.info.path, header = T, row.names = 1)

sig.count <- read.table(sig.count.path, header = T, row.names = 1)[, rownames(smp.info)] %>%
    t() %>%
    data.matrix()

if (!all(rownames(sig.count) == rownames(smp.info))) {
    stop("test.x and test.y not consistent")
}

print(paste0("AUC of ", trained.mdl.path, ": ", round(evalModel(trained.mdl, sig.count, smp.info$condition, model.name), 3)))
