rm(list = ls())

library(Matrix)
library(glmnet)
library(magrittr)
library(randomForest)

cmdArgs <- commandArgs(trailingOnly = T)
mdl.fit.path <- cmdArgs[1]
# mdl.fit.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/relapse/train0/lrc/fitted-model-ridge.rds"
new.mat.path <- cmdArgs[2]
# new.mat.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/relapse/test0/contig-counts-ridge.lrc.tsv"
out.path <- cmdArgs[3]
if (length(cmdArgs) == 4) {
    smp.info.path <- cmdArgs[4]
    # smp.info.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/relapse/sampleshuf.test0.tsv"
} else {
    smp.info.path <- NULL
}

write(mdl.fit.path, stderr())
write(new.mat.path, stderr())
write(out.path, stderr())
write(smp.info.path, stderr())

mdl.fit <- readRDS(mdl.fit.path)
new.mat <- read.table(new.mat.path, header = T, row.names = 1) %>% t()

features <- (coef(mdl.fit) %>% row.names())[-1]
mdl.name <- "glmnet"
if (is.null(features)) {
    features <- row.names(mdl.fit$importance)
    mdl.name <- "rf"
}
new.mat <- new.mat[, features]

pred.res <- data.frame("sample" = rownames(new.mat),
                       "prob.1" = predict(mdl.fit, new.mat, 
                                          type = ifelse(mdl.name == "glmnet", yes = "response", no = "prob")) %>% as.numeric(),
                       "pred.class" = predict(mdl.fit, new.mat, type = "class") %>% as.character())

write.table(pred.res, file = out.path, row.names = F, col.names = T, sep = "\t", quote = F)

if (!is.null(smp.info.path)) {
    smp.info <- read.table(smp.info.path, header = T, row.names = 1)
    smp.info <- smp.info[as.character(pred.res$sample), ]
    levels(pred.res$pred.class) <- levels(smp.info) # make the same levels between prediction and reality
    
    cat(paste0("balanced_accuracy\t", mdl.fit.path, "\t", 
               mlr3measures::bacc(truth = smp.info, response = pred.res$pred.class), "\n"))
    
    pred4auc <- ROCR::prediction(labels = smp.info, predictions = pred.res$prob.1)
    auc.pred <- ROCR::performance(pred4auc, measure = "auc")@y.values[[1]]
    cat(paste0("AUC\t", mdl.fit.path, "\t", auc.pred, "\n"))
}
