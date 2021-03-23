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
library(randomForest)
library(parallel)

cmdArgs <- commandArgs(trailingOnly = TRUE)
feature.mat.path <- cmdArgs[1]
# feature.mat.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/e_gene_level/b_rank_res/train0-lrc/top5000-gene-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/b_splitCV/sampleshuf.train0.tsv"
out.dir <- cmdArgs[3]
# out.dir <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0"
model.name <- cmdArgs[4]
# model.name <- "lasso"
if (length(cmdArgs) == 5) {
    nf.vect.path <- cmdArgs[5]
    # nf.vect.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/e_gene_level/b_rank_res/train0-lrc/smp-nf.tsv"
} else {
    nf.vect.path <- NA
}

print(feature.mat.path)
print(smp.info.path)
print(out.dir)
print(model.name)
print(nf.vect.path)

extractStableFeatures <- function(x, y, thres, num.runs, glmnet.alpha) {
    n <- nrow(x) 
    p <- ncol(x)
    features <- colnames(x) 
    cv.glm <- cv.glmnet(x = x, y = y, family = "binomial", alpha = glmnet.alpha, type.measure = "deviance")
    cl <- makeCluster(4)
    stabsel <- function(i) {
        cat("+")
        b_sort <- sort(sample(1 : n, round(3 * n/4)))
        out <- glmnet(x[b_sort, ], y[b_sort], family = "binomial",
                      lambda=cv.glm$lambda.1se, alpha = glmnet.alpha, standardize = FALSE)
        return(tabulate(which(out$beta[, ncol(out$beta)] != 0), p))
    }
    clusterEvalQ(cl, expr = c(library(glmnet)))
    clusterExport(cl, c('stabsel', 'num.runs', 'n', 'glmnet', 'cv.glm', 'x', 'y', 'p'),
                  envir = environment())
    res.cum <- Reduce("+", parLapply(cl, 1 : num.runs, stabsel))
    stopCluster(cl)
    prob.sel <- res.cum / num.runs
    feature.stabsel <- features[prob.sel >= thres]
    return(feature.stabsel)
}

fitModel <- function(train.x, train.y, mdl.name) {
    if (mdl.name == "lasso" || mdl.name == "elasticnet") {
        a <- ifelse(mdl.name == "lasso", yes = 1, no = 0.5) # if elastic-net, then alpha = 0.5
        feature.sel <- extractStableFeatures(x = train.x, y = train.y, thres = 0.5, num.runs = 2e3, glmnet.alpha = a)
        train.x <- train.x[, feature.sel]
        mdl.fit <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = a, type.measure = "deviance")
    } else if (mdl.name == "ridge") {
        mdl.fit <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = 0, type.measure = "deviance")
    } else if (mdl.name == "randomforest") {
        mdl.fit <- randomForest(x = train.x, y = train.y)
    } else {
        stop(paste0("unknown model name to fit: ", mdl.name, " (acceptable options: lasso, elasticnet, ridge, randomforest)"))
    }
    return(mdl.fit)
}

print(paste0("============> ", model.name, " <============"))

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

mdl.fit <- fitModel(train.x = feature.mat, train.y = smp.info$condition, mdl.name = model.name)
saveRDS(mdl.fit, file = paste0(out.dir, "/fitted-model-", model.name, ".rds"))

# Model Tuning Figure
png(paste0(out.dir, "/model-tuning-", model.name, ".png"), width = 800, height = 500)
par(mar = c(6, 5, 5, 3))
plot(mdl.fit, cex.axis = 2, cex.lab = 2, main = NULL)
title(paste0("model tuning: ", model.name), line = 3, cex.main = 2)
dev.off()

# Signatures and Model Assessment
if (model.name == "lasso" || model.name == "elasticnet" || model.name == "ridge") {
    sig.features <- rownames(coefficients(mdl.fit))[-1]
    if (length(sig.features) == 0) {
        print("no feature selected!")
    } else {
        auc.mdl <- assess.glmnet(mdl.fit, 
                                 newx = data.matrix(feature.mat[, sig.features]), newy = smp.info$condition,
                                 family = "binomial", s = "lambda.1se")$auc
        print(paste("AUC:", round(auc.mdl, 3)))
    }
} else {
    print(mdl.fit$confusion)
    sig.features <- rownames(importance(mdl.fit))
}

sig.features.fa <- cbind(paste0(">sig_feature_", seq(1 : length(sig.features))), 
                        matrix(sig.features, length(sig.features), byrow = T)) %>%
    t() %>%
    as.vector()
write.table(sig.features.fa, file = paste0(out.dir, "/signatures-", model.name, ".fa"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

print("")

