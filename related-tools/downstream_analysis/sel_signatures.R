###############################################
# A general GLMNET feature selection protocol #
# =========================================== #
# written by:                                 #
#     Thi Ngoc Ha Nguyen                      #
#     Haoliang Xue                            #
###############################################

rm(list = ls())

set.seed(91400)

library(stringr)
library(magrittr)
library(Matrix)
library(glmnet)
library(parallel)

cmdArgs <- commandArgs(trailingOnly = TRUE)
contig.mat.path <- cmdArgs[1]
# contig.mat.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/dekupl-run.train0/merged-diff-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.train0.tsv"
out.sig.path <- cmdArgs[3]
# out.sig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/signatures.fa"
out.mdl.path <- cmdArgs[4]
# out.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/model.rds"
out.fig.path <- cmdArgs[5]
# out.fig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/model.train0/lambda.tune.png"
glmnet.alpha <- cmdArgs[6] %>% as.numeric()
# glmnet.alpha <- 0
if (length(cmdArgs) == 7) {
    nf.vect.path <- cmdArgs[7]
    # nf.vect.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/c_kamrat_contigs/kamrat-snr.train0/smp-nf.tsv"
} else {
    nf.vect.path <- NA
}

extractSignatureStb <- function(train.x, train.y, thres, num.runs, a) {
    n <- nrow(train.x) 
    p <- ncol(train.x)
    features <- colnames(train.x) 
    cv.glm <- cv.glmnet(x = train.x, y = train.y, 
                        family = "binomial", alpha = a, type.measure = "deviance")
    cl <- makeCluster(4)
    stabsel <- function(i) {
        cat("+")
        b_sort <- sort(sample(1 : n, round(3 * n/4)))
        out <- glmnet(train.x[b_sort, ], train.y[b_sort], family = "binomial",
                      lambda=cv.glm$lambda.1se, alpha = a, standardize = FALSE)
        return(tabulate(which(out$beta[, ncol(out$beta)] != 0), p))
    }
    clusterEvalQ(cl, expr = c(library(glmnet)))
    clusterExport(cl, c('stabsel', 'num.runs', 'n', 'glmnet', 'cv.glm', 'train.x', 'train.y', 'p'),
                  envir = environment())
    res.cum <- Reduce("+", parLapply(cl, 1 : num.runs, stabsel))
    stopCluster(cl)
    prob.sel <- res.cum / num.runs
    feature.stabsel <- features[prob.sel >= thres]
    return(feature.stabsel)
}

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
    for (i in 1 : nrow(contig.mat)) {
        contig.mat[i, ] <- contig.mat[i, ] * nf.vect[rownames(contig.mat)[i], 1]
    }
}

if (!all(rownames(contig.mat) == smp.info$sample)) {
    stop("train.x and train.y not consistent")
}

contig.sel <- extractSignatureStb(train.x = contig.mat, train.y = smp.info$condition, 
                                  thres = 0.5, num.runs = 2e3, a = glmnet.alpha)

cv.mdl.postsel <- cv.glmnet(x = data.matrix(contig.mat[, contig.sel]), y = smp.info$condition,
                            family = "binomial", alpha = glmnet.alpha, type.measure = "deviance")
png(out.fig.path, width = 800, height = 500)
par(mar = c(6, 5, 5, 3))
plot(cv.mdl.postsel, cex.axis = 2, cex.lab = 2)
title(paste0("lambda tuning: alpha = ", glmnet.alpha), line = 3, cex.main = 2)
dev.off()
auc.postsel <- assess.glmnet(cv.mdl.postsel, 
                             newx = data.matrix(contig.mat[, contig.sel]), newy = smp.info$condition,
                             family = "binomial", s = "lambda.1se")$auc
print(paste("AUC:", round(auc.postsel, 3)))
saveRDS(cv.mdl.postsel, file = out.mdl.path)

sig.contigs.fa <- cbind(paste0(">sig_contig_", seq(1 : length(contig.sel))), 
                        matrix(contig.sel, length(contig.sel), byrow = T)) %>%
    t() %>%
    as.vector()
write.table(sig.contigs.fa, file = out.sig.path, quote = FALSE, row.names = FALSE, col.names = FALSE)
