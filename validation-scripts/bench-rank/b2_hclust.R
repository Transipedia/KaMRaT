rm (list = ls())

library(stringr)
library(ggdendro)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-rank/"
rank.dir <- paste0(work.dir, "/kamrat_res")

dlist <- c("ft20000-es10-pout0", "ft20000-es1.5-pout0", "ft20000-es1.5-pout20", "ft200000-es1.5-pout20")

for (d in dlist) {
    score.all <- NULL
    print(d)
    for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "bayes", "lr")) {
        f <- paste0(rank.dir, "/", d, "/ranked.", m, ".tsv")
        score.res <- read.table(f, header = T)[, c(1, 2)]
        score.res[, 2] <- 1 : nrow(score.res)
        colnames(score.res) <- c("feature", m)
        if (is.null(score.all)) {
            score.all <- score.res
        } else {
            score.all <- merge(score.all, score.res, by = "feature")
        }
    }
    cor.mat <- cor(score.all[, -1], method = "pearson")
    dist.mat <- as.dist((1 - cor.mat) / 2)
    p <- ggdendrogram(hclust(dist.mat, method = "complete"), 
                      rotate = T, theme_dendro = F, size = 1) +
        ylim(0, 0.5) +
        theme(text = element_text(size = 18, family = "Arial"),
              axis.text.x = element_text(hjust = 0.5),
              axis.title = element_blank(),
              plot.margin = margin(r = 0.1, t = 0.1, b = 0.1, l = 0.1, unit = "cm"))
    ggsave(filename = paste0(work.dir, "/", d, "-hclust.svg"), plot = p, width = 3.5, height = 2)
}
