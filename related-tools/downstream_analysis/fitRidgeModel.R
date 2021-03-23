########################################
# A general feature selection protocol #
# ==================================== #
# written by:                          #
#     Thi Ngoc Ha Nguyen               #
#     Haoliang Xue                     #
########################################

rm(list = ls())

set.seed(91400)

library(stringr)
library(magrittr)
library(Matrix)
library(glmnet)

cmdArgs <- commandArgs(trailingOnly = TRUE)
feature.mat.path <- cmdArgs[1]
# feature.mat.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/risk/train0/lrc/contig-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/risk/sampleshuf.train0.tsv"
out.dir <- cmdArgs[3]
# out.dir <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0"
ft.maxnum <- cmdArgs[4]
# ft.maxnum <- 100 %>% as.integer()
if (length(cmdArgs) == 5) {
    nf.vect.path <- cmdArgs[5]
    # nf.vect.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/risk/train0/lrc/smp-nf2.tsv"
} else {
    nf.vect.path <- NA
}

print(feature.mat.path)
print(smp.info.path)
print(out.dir)
print(ft.maxnum)
print(nf.vect.path)

fitModel <- function(train.x, train.y, feature_num) {
    return(mdl.fit)
}

smp.info <- read.table(smp.info.path, header = T)
smp.info <- smp.info[, c("sample", "condition")]

feature.mat <- read.table(feature.mat.path, header = T, row.names = 1)[, as.character(smp.info$sample)] %>%
    t() %>%
    data.matrix()

if (!is.na(nf.vect.path)) {
    print("    normalizing count table...")
    nf.vect <- read.table(nf.vect.path, header = FALSE, row.names = 1)
    for (s in rownames(nf.vect)) {
        feature.mat[s, ] <- feature.mat[s, ] * nf.vect[s, 1]
    }
}

if (!all(rownames(feature.mat) == smp.info$sample)) {
    stop("train.x and train.y not consistent")
}

mdl.fit <- cv.glmnet(x = feature.mat, y = smp.info$condition, family = "binomial", alpha = 0, type.measure = "deviance")
ft.ipt <- abs(coefficients(mdl.fit)[-1, ]) %>% as.data.frame()
top.ft.mat <- feature.mat[, rownames(ft.ipt)[order(ft.ipt$., decreasing = T)[1 : ft.maxnum]]]
mdl.fit <- cv.glmnet(x = top.ft.mat, y = smp.info$condition, family = "binomial", alpha = 0, type.measure = "deviance")

saveRDS(mdl.fit, file = paste0(out.dir, "/fitted-model-", ft.maxnum, ".rds"))

# Model Tuning Figure
png(paste0(out.dir, "/model-tuning-", ft.maxnum, ".png"), width = 800, height = 500)
par(mar = c(6, 5, 5, 3))
plot(mdl.fit, cex.axis = 2, cex.lab = 2, main = NULL)
title(paste0("model tuning: ", ft.maxnum), line = 3, cex.main = 2)
dev.off()

# Signatures and Model Assessment
sig.features <- rownames(coefficients(mdl.fit))[-1]
if (length(sig.features) == 0) {
    print("no feature selected!")
} else {
    auc.mdl <- assess.glmnet(mdl.fit, 
                             newx = data.matrix(feature.mat[, sig.features]), newy = smp.info$condition,
                             family = "binomial", s = "lambda.1se")$auc
    print(paste("AUC:", round(auc.mdl, 3)))
}

sig.features.fa <- cbind(paste0(">sig_feature_", seq(1 : length(sig.features))), 
                        matrix(sig.features, length(sig.features), byrow = T)) %>%
    t() %>%
    as.vector()
write.table(sig.features.fa, file = paste0(out.dir, "/signatures-", ft.maxnum, ".fa"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print("")

