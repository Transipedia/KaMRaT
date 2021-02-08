rm(list = ls())

library(stringr)
library(magrittr)
library(ModelMetrics) # F1 score
library(e1071) # Naive Bayes, SVM
library(SVMMaj) # Hinge
library(FSelectorRcpp) # Information Gain

cmdArgs <- commandArgs(trailingOnly = T)
count.path <- cmdArgs[1]
# count.path <- "/home/haoliang.xue/development/kamrat-test/data/test-counts.tsv"
smp.info.path <- cmdArgs[2]
# smp.info.path <- "/home/haoliang.xue/development/kamrat-test/data/test-sample-conditions.tsv"
out.dir <- cmdArgs[3]
# out.dir <- ""
no.norm <- cmdArgs[4] %>% as.logical()
# no.norm <- FALSE

calcRelatSD <- function(count_tab, smp_nf) {
    row.mean <- apply(count_tab[, names(smp_nf)], MARGIN = 1, function(r) mean(r * smp_nf))
    row.sd <- apply(count_tab[, names(smp_nf)], MARGIN = 1, function(r) sd(r * smp_nf))
    row.rsd <- lapply(1 : nrow(count_tab), 
                      FUN = function(i) ifelse(row.mean[i] > 1, yes = row.sd[i] / row.mean[i], no = row.sd[i])) %>%
        unlist()
    return(row.rsd)
}

calcTtestPvalue <- function(count_tab, smp_info, smp_nf) {
    cond.list <- unique(smp_info[, 1]) %>% 
        as.character()
    cond1.cols <- names(count_tab)[smp_info[names(count_tab), 1] == cond.list[1]]
    cond2.cols <- names(count_tab)[smp_info[names(count_tab), 1] == cond.list[2]]
    pval.raw <- apply(count_tab, MARGIN = 1, 
                      function(r) t.test(x = log(r[cond1.cols] * smp_nf[cond1.cols] + 1), 
                                         y = log(r[cond2.cols] * smp_nf[cond2.cols] + 1),
                                         alternative = "two.sided")$p.value)
    pval.adj <- p.adjust(pval.raw, method = "BH")
    return(pval.adj)
}

calcSNR <- function(count_tab, smp_info, smp_nf) {
    cond.list <- unique(smp_info[, 1]) %>% 
        as.character()
    cond1.cols <- names(count_tab)[smp_info[names(count_tab), 1] == cond.list[1]]
    cond2.cols <- names(count_tab)[smp_info[names(count_tab), 1] == cond.list[2]]
    row.meanA <- apply(count_tab[, cond1.cols], MARGIN = 1, function(r) mean(r * smp_nf[cond1.cols]))
    row.meanB <- apply(count_tab[, cond2.cols], MARGIN = 1, function(r) mean(r * smp_nf[cond2.cols]))
    row.sdA <- apply(count_tab[, cond1.cols], MARGIN = 1, function(r) sd(r * smp_nf[cond1.cols]))
    row.sdB <- apply(count_tab[, cond2.cols], MARGIN = 1, function(r) sd(r * smp_nf[cond2.cols]))
    row.snr <- (row.meanA - row.meanB) / (row.sdA + row.sdB)
    return(row.snr)
}

calcLRCF1 <- function(count_tab, smp_info, smp_nf) {
    for (s in names(smp_nf)) {
        count_tab[, s] <- count_tab[, s] * smp_nf[s]
    }
    smp_info$label <- as.numeric(smp_info$V2) - 1
    tmp.df <- merge(x = t(count_tab), y = smp_info, by = "row.names")
    row.f1 <- lapply(rownames(count_tab),
                     function(r) glm(formula = tmp.df$label ~ tmp.df[, r],
                                     family = binomial(link = "logit"))$fitted.values %>%
                         f1Score(actual = tmp.df$label, cutoff = 0.5)) %>%
        unlist()
    return(row.f1)
}

calcNBCF1 <- function(count_tab, smp_info, smp_nf) {
    for (s in names(smp_nf)) {
        count_tab[, s] <- count_tab[, s] * smp_nf[s]
    }
    smp_info$label <- as.numeric(smp_info$V2) - 1
    tmp.df <- merge(x = t(count_tab), y = smp_info, by = "row.names")
    row.f1 <- lapply(rownames(count_tab),
                     function(r) naiveBayes(x = as.matrix(tmp.df[, r]), y = factor(tmp.df$label)) %>%
                         predict(newdata = as.matrix(tmp.df[, r]), type = "class") %>%
                         f1Score(actual = tmp.df$label, cutoff = 1)) %>%
        unlist()
    return(row.f1)
}

calcSVMHingeLoss <- function(count_tab, smp_info, smp_nf) {
    hingefunc <- getHinge(hinge = "absolute", delta = 1)
    for (s in names(smp_nf)) {
        count_tab[, s] <- count_tab[, s] * smp_nf[s]
        count_tab[, s] <- (count_tab[, s] - mean(count_tab[, s])) / sd(count_tab[, s])
    }
    smp_info$label <- ifelse(as.numeric(smp_info$V2) == 1, yes = -1, no = 1)
    tmp.df <- merge(x = t(count_tab), y = smp_info, by = "row.names")
    hingeloss.list <- lapply(rownames(count_tab),
                             function(r) svm(x = as.matrix(tmp.df[, r]), y = tmp.df$label, 
                                             kernel = "linear", cost = 1) %>%
                                 predict(newdata = as.matrix(tmp.df[, r])) %>%
                                 hingefunc(y = tmp.df$label))
    row.hingeloss <- lapply(hingeloss.list,  function(r) sum(r$loss)) %>%
        unlist()
    return(row.hingeloss)
}

calcInfoGain <- function(count_tab, smp_info, smp_nf) {
    for (s in names(smp_nf)) {
        count_tab[, s] <- count_tab[, s] * smp_nf[s]
    }
    row.infogain <- apply(count_tab, MARGIN = 1,
                          function(r) information_gain(x = as.data.frame(r), 
                                                       y = as.vector(smp.info$V2), 
                                                       type = "infogain")$importance) %>%
        unlist()
    return(row.infogain)
}

smp.info <- read.table(smp.info.path, header = F, row.names = 1)
count.tab <- read.table(count.path, header = T, row.names = 1)[, rownames(smp.info)]

if (!no.norm) {
    smp.sum <- colSums(count.tab)
    smp.nf <- mean(smp.sum) / smp.sum
} else {
    smp.nf <- rep(1, nrow(smp.info))
    names(smp.nf) <- row.names(smp.info)
}

feature.score <- data.frame(row.names = rownames(count.tab), 
                            "relat_sd" = calcRelatSD(count.tab, smp.nf),
                            "ttest.padj" = calcTtestPvalue(count.tab, smp.info, smp.nf),
                            "snr" = calcSNR(count.tab, smp.info, smp.nf),
                            "lrc.f1" = calcLRCF1(count.tab, smp.info, smp.nf),
                            "nbc.f1" = calcNBCF1(count.tab, smp.info, smp.nf),
                            "svm.hingeloss" = calcSVMHingeLoss(count.tab, smp.info, smp.nf),
                            "info.gain" = calcInfoGain(count.tab, smp.info, smp.nf))

write.table(feature.score, file = paste0(out.dir, "/feature-scores.byR.", ifelse(no.norm, yes = "raw", no = "norm"), ".tsv"), 
            col.names = T, row.names = T, quote = F, sep = "\t")
