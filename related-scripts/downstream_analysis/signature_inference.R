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
library(caret)

cmdArgs <- commandArgs(trailingOnly = TRUE)
top.contig.path <- cmdArgs[1]
top.contig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/data/merged-diff-counts.tsv"
sample.info.path <- cmdArgs[2]
sample.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/data/sample-info.txt"
out.sig.path <- cmdArgs[3]
out.sig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/signatures.fa"
out.mdl.path <- cmdArgs[4]
out.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/model.rds"
if (length(cmdArgs) == 5) {
    nf.vect.path <- cmdArgs[5]
    # nf.vect.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-test/sample-nf-snr.tsv"
} else {
    nf.vect.path <- NA
}

ExtractSignatureStb <- function(df, thres, num.runs){
    x.train <- df[, !(names(df) %in% c("sample", "condition"))] %>% data.matrix()
    y.train <- df$condition
    feature.names <- colnames(x.train)
    model.cv.glm <- cv.glmnet(x = x.train, y = y.train, 
                              family = "binomial", alpha = 1, type.measure = "mse")
    StabSel <- function(i) {
        n <- nrow(x.train)
        b_sort <- sort(sample(1 : n, round(3 * n / 4)))
        out <- glmnet(x.train[b_sort, ], y.train[b_sort],
                      lambda = model.cv.glm$lambda.1se, family = "binomial", 
                      alpha = 1, standardize = FALSE)
        return(as.numeric(out$beta[, ncol(out$beta)] != 0))
    }
    res.cum <- lapply(1 : num.runs, FUN = StabSel) %>%
        do.call(what = rbind) %>% 
        colSums()
    prob.sel <- res.cum / num.runs
    # plot(sort(prob.sel))
    feature.stabsel <- feature.names[prob.sel >= thres]
    return (feature.stabsel)
}

print(paste0("Top contig path:", top.contig.path))

sample_info <- read.table(sample.info.path, header = F)
names(sample_info) <- c("sample", "condition")

top_contigs <- read.table(top.contig.path, header = T, sep = "\t")
rownames(top_contigs) <- top_contigs$contig
top_contigs <- top_contigs[, names(top_contigs) %in% sample_info$sample]

if (!is.na(nf.vect.path)) {
    print("    normalizing count table...")
    nf.vect <- read.table(nf.vect.path, header = FALSE, row.names = 1)
    for (i in 1 : nrow(nf.vect)) {
        top_contigs[, row.names(nf.vect)[i]] <- top_contigs[, row.names(nf.vect)[i]] * nf.vect[i, 1]
    }
}

feature.tab <- as.data.frame(t(top_contigs))
feature.tab$sample <- row.names(feature.tab)

data.tab <- merge(x = feature.tab, y = sample_info, by = "sample", all.x = TRUE)

thres.lasso = 0.5 # threshold for Lasso logistic regression
sig.contigs <- ExtractSignatureStb(data.tab, thres = thres.lasso, num.runs = 2e3)

# Train model with selected signatures
data.tab.sel <- data.tab[, c(sig.contigs, "condition")]
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 20,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = TRUE)
mdl.postsel <- train(condition ~ ., data = data.tab.sel,
                     method = "glm",
                     metric = "ROC",
                     family = "binomial",
                     preProc = c("center", "scale"),
                     trControl = ctrl)  

# Save model and contig fasta
saveRDS(mdl.postsel, file = out.mdl.path)

contig.names <- paste0(">contig_", seq(1 : length(sig.contigs)))
sig.contigs.fa <- as.vector(t(cbind(contig.names, matrix(sig.contigs, length(sig.contigs), byrow = T))))
out.filename <- paste0(out.sig.path) %>%
    str_replace(pattern = ".tsv", replacement = ".fa")
write.table(sig.contigs.fa, file = out.filename, 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
