rm (list = ls())

library(stringr)
library(magrittr)
library(parallel)
library(foreach)
library(tidyr)

work.dir <- "../"

tab.in <- read.table(gzfile(paste0(work.dir, "/data/kmer-counts.subset4toy.tsv.gz")), 
                     header = T, row.names = 1)
for (i in colnames(tab.in)) {
    tab.in[, i] <- tab.in[, i] / sum(tab.in[, i]) * 10000000
}

evalRow <- function(X) {
    library(magrittr)
    library(MLmetrics) # Accuracy
    library(e1071) # NBC, LR
    library(DescTools) # Entropy
    
    work.dir <- "../"
    
    smp.condi.categ <- read.table(paste0(work.dir, "/data/sample-condition.toy.tsv"))
    names(smp.condi.categ) <- c("sample", "condition")
    smp.condi.categ$condition <- as.numeric(smp.condi.categ$condition == unique(smp.condi.categ$condition)[2])
    
    smp.condi.cntnu <- read.table(paste0(work.dir, "/data/sample-condition.toy2.tsv"))
    names(smp.condi.cntnu) <- c("sample", "condition")
    smp.condi.cntnu$condition <- as.numeric(smp.condi.cntnu$condition)
    
    X.df <- data.frame("sample" = names(X), "count" = as.numeric(X))
    s <- sd(X.df$count) # sd
    rsd1 <- s / max(c(1, mean(X.df$count))) # sd/mean
    rsd2 <- s / max(c(1, min(X.df$count))) # sd/min
    etp <- Entropy(X.df$count + 1) # entropy
    
    df.categ <- merge(x = X.df, y = smp.condi.categ, by = "sample")
    x1 <- df.categ$count[df.categ$condition == 0]
    x2 <- df.categ$count[df.categ$condition == 1]
    
    if (s != 0) {
        ttest.praw <- t.test(x = log2(x1 + 1), y = log2(x2 + 1), 
                             alternative = "two.sided")$p.value # to calculate ttest.padj
        log2fc <- abs(mean(log2(x1 + 1)) - mean(log2(x2 + 1)))
        ttest.pi <- -log10(ttest.praw) * log2fc # ttest.pi
        snr <- (mean(x1) - mean(x2))/(sd(x1) + sd(x2)) # SNR
    } else {
        ttest.praw <- 1
        ttest.pi <- 0
        snr <- 0
    }
    
    dids <- lapply(unique(df.categ$condition), 
                   FUN = function(i) {
                       grp.max <- max(df.categ$count[df.categ$condition == i])
                       grp.dist <- df.categ$count[df.categ$condition != i] - grp.max
                       grp.dist[grp.dist < 0] <- 0
                       grp.dids <- sum(sqrt(grp.dist))
                   }) %>% unlist() %>% max()
    
    df.categ$count <- (df.categ$count - mean(df.categ$count)) / sd(df.categ$count) # standardization for ML
    
    nbc.acc <- naiveBayes(x = as.matrix(df.categ$count), y = factor(df.categ$condition)) %>% 
        predict(newdata = as.matrix(df.categ$count), type = "class") %>%
        Accuracy(y_true = df.categ$condition) # Bayes
    
    if (s > 0) {
        lr.pred <- glm(formula = condition ~ count, data = df.categ, family = binomial) %>%
            predict(newdata = df.categ, type = "response")
        lr.acc <- Accuracy(y_pred = factor(ifelse(lr.pred < 0.5, yes = 0, no = 1), levels = c(0, 1)), 
                           y_true = df.categ$condition) # LR
    } else {
        lr.acc <- 0
    }
    
    df.cntnu <- merge(x = X.df, y = smp.condi.cntnu, by = "sample")
    
    ps.cor <- ifelse(s == 0, 
                     no = cor(df.cntnu$count, df.cntnu$condition, method = "pearson"),
                     yes = 0)
    
    sp.cor <- ifelse(s == 0,
                     no = cor(df.cntnu$count, df.cntnu$condition, method = "spearman"),
                     yes = 0)
    
    return(data.frame("ttest.praw" = ttest.praw,
                      "ttest.pi" = ttest.pi,
                      "SNR" = snr,
                      "LR.acc" = lr.acc,
                      "DIDS.score" = dids,
                      "Bayes.acc" = nbc.acc,
                      "pearson" = ps.cor,
                      "spearman" = sp.cor,
                      "sd" = s,
                      "rsd1" = rsd1,
                      "rsd2" = rsd2,
                      "entropy" = etp))
}

cl <- makeCluster(4)
eval.res <- parApply(cl = cl, tab.in, MARGIN = 1, FUN = evalRow) %>%
    do.call(what = rbind)
stopCluster(cl)
eval.res$ttest.padj <- p.adjust(eval.res$ttest.praw, method = "BH") # ttest.padj

eval.res$tag <- rownames(eval.res)
eval.res.long <- pivot_longer(eval.res, cols = -tag, names_to = "method", values_to = "score.R")

kamrat.res <- NULL
for (f in dir(paste0(work.dir, "/output/kamrat-rank/"))) {
    if (is.null(kamrat.res)) {
        kamrat.res <- read.table(paste0(work.dir, "/output/kamrat-rank/", f), header = T)[, 1 : 2]
    } else {
        kamrat.res <- merge(x = kamrat.res,
                            y = read.table(paste0(work.dir, "/output/kamrat-rank/", f), 
                                       header = T)[, 1 : 2],
                            by = "tag")
    }
}

kamrat.res.long <- pivot_longer(kamrat.res, cols = -tag, names_to = "method", values_to = "score.kamrat")

cmp.res <- merge(kamrat.res.long, eval.res.long, by = c("tag", "method"))
cmp.res$diff <- abs(cmp.res$score.kamrat - cmp.res$score.R)
cmp.res$relat.diff <- abs(cmp.res$diff / cmp.res$score.R)
