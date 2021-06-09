rm(list = ls())

library(stringr)
library(magrittr)
library(MLmetrics) # Accuracy
library(e1071) # Naive Bayes, SVM
library(SVMMaj) # Hinge
library(FSelectorRcpp) # Information Gain

count.path <- "/home/haoliang.xue/development/new-kamrat-test/data/luad-diff-counts-up-samples.tsv"
smp.info.path <- "/home/haoliang.xue/development/new-kamrat-test/data/luad-sample-conditions.tsv"
rank.dir <- "/home/haoliang.xue/development/new-kamrat-test/rank-res/"
out.path <- "/home/haoliang.xue/development/new-kamrat-test/rank-res/feature-scores.byR.tsv"

smp.info <- read.table(smp.info.path, header = F, row.names = 1)
count.tab <- read.table(count.path, header = T, row.names = 1)[, rownames(smp.info)]
smp.sum <- colSums(count.tab)
for (x in names(smp.sum)) {
    count.tab[, x] <- count.tab[, x] / smp.sum[x] * 1000000000
}

evalRow <- function(X, y) {
    d <- merge(data.frame("count" = X, "row.names" = names(X)), y, by = "row.names")
    names(d) <- c("sample", "count", "condition")
    d$label <- as.numeric(d$condition) - 1
    # t-test non-adjusted p-value
    row.praw <- t.test(x = log(d$count[d$label == 0] + 1), y = log(d$count[d$label == 1] + 1), 
                       alternative = "two.sided")$p.value
    # SNR
    row.snr <- (mean(d$count[d$label == 0]) - mean(d$count[d$label == 1]))/(sd(d$count[d$label == 0]) + sd(d$count[d$label == 1]))
    # F1 score and accuracy of naive Bayes classifier
    row.nbc.acc <- naiveBayes(x = as.matrix(d$count), y = factor(d$label)) %>% 
        predict(newdata = as.matrix(d$count), type = "class") %>%
        Accuracy(y_true = d$label)
    # F1 score and accuracy of logistic regression
    lr.pred <- glm(formula = label ~ count, data = d, family = binomial) %>%
        predict(newdata = d, type = "response")
    row.lr.acc <- Accuracy(y_pred = factor(ifelse(lr.pred < 0.5, yes = 0, no = 1), levels = c(0, 1)), 
                           y_true = d$label)
    row.svm.acc <- svm(x = as.matrix(d$count), y = factor(d$label), kernel = "linear", cost = 1) %>%
        predict(newdata = as.matrix(d$count)) %>%
        Accuracy(y_true = d$label)
    return(data.frame("ttest.praw" = row.praw,
                      "snr" = row.snr,
                      "nbc.acc" = row.nbc.acc,
                      "lr.acc" = row.lr.acc,
                      "svm.acc" = row.svm.acc))
}

feature.score <- apply(count.tab, MARGIN = 1, FUN = function(r) evalRow(r, smp.info))
write.table(feature.score, file = out.path, 
            col.names = T, row.names = T, quote = F, sep = "\t")
