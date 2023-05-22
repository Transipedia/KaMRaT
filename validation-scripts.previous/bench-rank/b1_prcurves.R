rm(list = ls())

library(stringr)
library(magrittr)
library(PRROC)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-rank/"
data.dir <- paste0(work.dir, "compcodeR_res/")
rank.dir <- paste0(work.dir, "kamrat_res/")

dlist <- c("ft20000-es10-pout0", "ft20000-es1.5-pout0", "ft20000-es1.5-pout20", "ft200000-es1.5-pout20")
ypos.list <- c(0.1, 0.2, 0.3, 0.4)

label.real <- lapply(dlist, 
                     FUN = function(d) {
                         df <- read.table(paste0(data.dir, d, "/varinfo.tsv"), 
                                          header = T)[, c("feature", "differential.expression")]
                         df$dataset <- d
                         return(df)
                         }) %>%
    do.call(what = rbind)

auc.all <- NULL
pr.all <- NULL
for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "bayes", "lr")) {
    print(m)
    for (i in 1 : 4) {
        d <- dlist[i]
        f <- paste0(rank.dir, "/", d, "/ranked.", m, ".tsv")
        score.res <- read.table(f, header = T)[, c(1, 2)]
        names(score.res) <- c("feature", "score")
        if (m == "snr") {
            score.res$score <- abs(score.res$score)
        } else if (m == "ttest.padj") {
            score.res$score <- -score.res$score
        }
        score.res$dataset <- d
        cmp.res <- merge(x = label.real, y = score.res, by = c("feature", "dataset"))
        pr <- pr.curve(scores.class0 = cmp.res$score, 
                       weights.class0 = cmp.res$differential.expression, curve = T, minStepSize = 0.05)
        pr.values <- as.data.frame(pr$curve[, -3])
        names(pr.values) <- c("recall", "precision")
        pr.values$rk.mthd <- m
        pr.values$dataset <- d
        pr.all <- rbind(pr.all, pr.values)
        auc.all <- rbind(auc.all, data.frame("auc" = pr$auc.integral, "rk.mthd" = m, "dataset" = d, "ypos" = ypos.list[i]))
    }    
}
write.table(auc.all, paste0(work.dir, "/PRauc-compcodeR.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
p <- ggplot() +
    geom_line(data = pr.all, aes(x = recall, y = precision, color = dataset), size = 1) +
    facet_wrap(rk.mthd ~ ., nrow = 2) +
    scale_color_manual(values = c("#b2abd2", "#fdb863", "#5e3c99", "#e66101")) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2))) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2))) +
    coord_fixed() +
    xlab("recall") +
    ylab("precision") +
    theme(text = element_text(size = 35, family = "Arial"), legend.position = "none")
ggsave(filename = paste0(work.dir, "/PRcurves-compcodeR.svg"), plot = p, width = 16, height = 12)
