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
library(glmnet)
library(parallel)
library(caret)

ExtractSignatureStb <- function(df, thres){
  X <- df[, !(names(df) %in% c("sample", "condition"))] %>% data.matrix()
  Y <- df$condition
  feature.names <- colnames(X) 
  model.cv.glm <- cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1, type.measure = "mse")
  cores <- detectCores()
  cl <- makeCluster(cores - 2)
  NUM_RUNS <- 2e3
  n <- nrow(X)
  p <- ncol(X)
  stabsel <- function(i) {
    cat("+")
    b_sort <- sort(sample(1 : n, round(3 * n / 4)))
    out <- glmnet(X[b_sort, ], Y[b_sort],
                  lambda = model.cv.glm$lambda.1se, family = "binomial", 
                  alpha = 1, standardize = FALSE)
    return(tabulate(which(out$beta[, ncol(out$beta)] != 0), p))
  }
  clusterEvalQ(cl, expr = c(library(glmnet)))
  clusterExport(cl, c('stabsel', 'frame', 'NUM_RUNS', 'n', 'glmnet', 'model.cv.glm', 'X', 'Y', 'p'), 
                envir = environment())
  res.cum <- Reduce("+", parLapply(cl, 1 : NUM_RUNS, stabsel))
  stopCluster(cl)
  prob.sel <- res.cum / NUM_RUNS
  # plot(sort(prob.sel))
  feature.stabsel <- feature.names[prob.sel >= thres]
  return (feature.stabsel)
}

cmdArgs <- commandArgs(trailingOnly = TRUE)
top.contig.path <- cmdArgs[1]
# top.contig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/merged-masked-count.snr20000.tsv"
sample.info.path <- cmdArgs[2]
# sample.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/data/sample-info.txt"
out.dir <- cmdArgs[3]
# out.dir <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results"
nf.vect.path <- cmdArgs[4]
# nf.vect.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-test/sample-nf-snr.tsv"

print(paste0("Top contig path:", top.contig.path))

out.filename <- paste0(out.dir, "/signature-", basename(top.contig.path)) %>%
  str_replace(pattern = ".tsv", replacement = ".fa")

sample_info <- read.table(sample.info.path, header = F)
names(sample_info) <- c("sample", "condition")

NUM_RUNS=100 # Number of time for sampling

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
sig.contigs <- ExtractSignatureStb(data.tab, thres = thres.lasso)

# Train model with selected signatures
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 20,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = TRUE)
glm.fit <- train(condition ~ ., data = data.tab[, c(sig.contigs, "condition")],
                method = "glm",
                metric = "ROC",
                family = "binomial",
                preProc = c("center", "scale"),
                trControl = ctrl)  
auc.mean <- round(mean(glm.fit$resample$ROC), 3)
print(paste("    mean AUC:", auc.mean))

# Create fa file
contig.names <- paste0(">contig_", seq(1 : length(sig.contigs)))
sig.contigs.fa <- as.vector(t(cbind(contig.names, matrix(sig.contigs, length(sig.contigs), byrow = T))))
  
write.table(sig.contigs.fa, file = out.filename, 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
