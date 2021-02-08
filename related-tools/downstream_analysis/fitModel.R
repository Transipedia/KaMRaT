rm(list = ls())

set.seed(91400)

library(stringr)
library(magrittr)
library(Matrix)
library(glmnet)
library(randomForest)

source("/store/USERS/haoliang.xue/development/KaMRaT/related-tools/downstream_analysis/selectFeatures.func.R")

cmdArgs <- commandArgs(trailingOnly = TRUE)
contig.mat.path <- cmdArgs[1]
# contig.mat.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/dekupl-run.train0/merged-diff-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.train0.tsv"
out.dir <- cmdArgs[3]
# out.dir <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0"
model.name <- cmdArgs[4]
# model.name <- "lasso"
if (length(cmdArgs) == 5) {
    nf.vect.path <- cmdArgs[5]
    # nf.vect.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/c_kamrat_contigs/kamrat-snr.train0/smp-nf.tsv"
} else {
    nf.vect.path <- NA
}

print(paste0("============> ", model.name, " <============"))

smp.info <- read.table(smp.info.path, header = T)
smp.info <- smp.info[, c("sample", "condition")]

contig.mat <- read.table(contig.mat.path, header = T)
rownames(contig.mat) <- contig.mat$contig
contig.mat <- contig.mat[, as.character(smp.info$sample)] %>%
    t() %>%
    data.matrix()

if (!is.na(nf.vect.path)) {
    print("    normalizing count table...")
    nf.vect <- read.table(nf.vect.path, header = FALSE, row.names = 1)
    for (s in rownames(nf.vect)) {
        contig.mat[s, ] <- contig.mat[s, ] * nf.vect[s, 1]
    }
}

if (!all(rownames(contig.mat) == smp.info$sample)) {
    stop("train.x and train.y not consistent")
}

mdl.fit <- fitModel(train.x = contig.mat, train.y = smp.info$condition, mdl.name = model.name)
saveRDS(mdl.fit, file = paste0(out.dir, "/fitted-model-", model.name, ".rds"))

# Model Tuning Figure
png(paste0(out.dir, "/model-tuning-", model.name, ".png"), width = 800, height = 500)
par(mar = c(6, 5, 5, 3))
plot(mdl.fit, cex.axis = 2, cex.lab = 2, main = NULL)
title(paste0("model tuning: ", model.name), line = 3, cex.main = 2)
dev.off()

# Signatures and Model Assessment
if (model.name == "lasso" || model.name == "elasticnet" || model.name == "ridge") {
    sig.contigs <- rownames(coefficients(mdl.fit))[-1]
    auc.mdl <- assess.glmnet(mdl.fit, 
                             newx = data.matrix(contig.mat[, sig.contigs]), newy = smp.info$condition,
                             family = "binomial", s = "lambda.1se")$auc
    print(paste("AUC:", round(auc.mdl, 3)))
} else {
    print(mdl.fit$confusion)
    sig.contigs <- rownames(importance(mdl.fit))
}

sig.contigs.fa <- cbind(paste0(">sig_contig_", seq(1 : length(sig.contigs))), 
                        matrix(sig.contigs, length(sig.contigs), byrow = T)) %>%
    t() %>%
    as.vector()
write.table(sig.contigs.fa, file = paste0(out.dir, "/signatures-", model.name, ".fa"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print("")