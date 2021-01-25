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

cmdArgs <- commandArgs(trailingOnly = TRUE)
top.contig.path <- cmdArgs[1]
# top.contig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/c_kamrat_contigs/kamrat-snr.train0/merged-ranked-masked-counts.top10000.tsv"
sample.info.path <- cmdArgs[2]
# sample.info.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/a_splitted_dataset/sampleshuf.train0.tsv"
out.sig.path <- cmdArgs[3]
# out.sig.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/signatures.fa"
out.mdl.path <- cmdArgs[4]
# out.mdl.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/new.results/model.rds"
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
                              family = "binomial", alpha = 1, type.measure = "mse") # mse ? deviance ?
    StabSel <- function(i) {
        n <- nrow(x.train)
        b_sort <- sort(sample(1 : n, round(3 * n / 4)))
        out <- glmnet(x = x.train[b_sort, ], y = y.train[b_sort], lambda = model.cv.glm$lambda.1se, 
                      family = "binomial", alpha = 1, standardize = FALSE) # standardize = FALSE ? TRUE ?
        return(as.numeric(out$beta[, ncol(out$beta)] != 0))
    }
    res.cum <- lapply(1 : num.runs, FUN = StabSel) %>%
        do.call(what = rbind) %>% 
        colSums()
    prob.sel <- res.cum / num.runs
    feature.stabsel <- feature.names[prob.sel >= thres]
    return (list("lambda.1se" = model.cv.glm$lambda.1se, "feature.sel" = feature.stabsel))
}

print(paste0("Top contig path:", top.contig.path))

sample_info <- read.table(sample.info.path, header = T)
sample_info <- sample_info[, c("sample", "condition")]

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
sel.res <- ExtractSignatureStb(data.tab, thres = thres.lasso, num.runs = 2e3)
lambda.1se <- sel.res[["lambda.1se"]]
contig.sel <- sel.res[["feature.sel"]]

# Train model with selected signatures
data.tab.sel <- data.tab[, c(contig.sel, "condition")]
mdl.postsel <- glmnet(x = data.matrix(data.tab[, contig.sel]),
                      y = data.tab$condition, type.measure = "auc", # mse ? auc ? deviance ?
                      lambda = lambda.1se,
                      family = "binomial", alpha = 1, standardize = TRUE)

# Save model and contig fasta
saveRDS(list("trained.mdl" = mdl.postsel, "lambda.1se" = lambda.1se), file = out.mdl.path)

contig.names <- paste0(">sig_contig_", seq(1 : length(contig.sel)))
sig.contigs.fa <- as.vector(t(cbind(contig.names, matrix(contig.sel, length(contig.sel), byrow = T))))
out.filename <- paste0(out.sig.path) %>%
    str_replace(pattern = ".tsv", replacement = ".fa")
write.table(sig.contigs.fa, file = out.filename, 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
