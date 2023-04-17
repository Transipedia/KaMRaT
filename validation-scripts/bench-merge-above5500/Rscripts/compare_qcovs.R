rm(list = ls())

#library(ggpubr)
library(ggplot2)

work.dir <- "../../../../revision/bench-merge-above5500/"
# work.dir <- "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"

read_align <- function(depth, merge_mode) {
    align.res <- read.table(paste0(work.dir, "depth_", depth, 
                                   "/ctg-aligned.",
                                   merge_mode, ".tsv"),
                            sep="\t", header = TRUE, stringsAsFactors=F)
    align.res <- align.res[align.res$qlen > 61, ]
    align.res$source <- merge_mode
    return(align.res)
}
align <- NULL
for (dpt in c(0.05, 0.2)) {
    for (m in c("none", "pearson_0.2", "spearman_0.2", "mac_0.2")) {
        cat(dpt, m, "\n")
        align.x <- read_align(dpt, m)
        align.x$depth <- dpt
        align.x$method <- m
        align <- rbind(align, align.x)
    }    
}


# pval.df <- NULL
# for(merge_mode1 in c("none", "pearson", "spearman", "mac")) {
#     for (merge_mode2 in c("none", "pearson", "spearman", "mac")) {
#         pval <- t.test(unlist(align[[merge_mode1]]["qcovs"]),
#                             unlist(align[[merge_mode2]]["qcovs"]),
#                             alternative = "two.sided")$p.value
#         pval.df <- rbind(pval.df,
#                          data.frame("merge.mode.1" = merge_mode1,
#                                     "merge.mode.2" = merge_mode2,
#                                     "pvalue.nominal" = pval))
#     }
# }
# pval.df$pvalue.adj <- p.adjust(pval.df$pvalue.nominal, method = "bonferroni")
# pval.df$significant <- pval.df$pvalue.adj < 0.05
# pval.df <- pval.df[order(pval.df$pvalue.adj, pval.df$pvalue.nominal, decreasing = FALSE), ]
# write.csv(pval.df, file = paste0(work.dir, "kamrat_res/merge-on-real/PRADtcga/qcovs_wx-test.csv"),
#           quote = FALSE, row.names = FALSE)

# p <- ggboxplot(rbind(align[["none"]], align[["mac"]], align[["pearson"]], align[["spearman"]]),
#                x = "source", y = "qcovs", color = "source",
#                add = "jitter", size = 1, notch = TRUE) + 

p <- ggplot(data = align, aes(x = paste(as.character(depth), method),
                         y = qcovs, fill = method)) +
    # geom_boxplot() +
    geom_jitter(colour = "black", size = 0.2, alpha = 0.05) +
    scale_fill_manual(values = c("none" = "#e66101",
                                 "mac_0.2" = "#fdb863",
                                 "pearson_0.2" = "#b2abd2",
                                 "spearman_0.2" = "#5e3c99")) +
    theme_bw() +
    theme(text = element_text(size = 15))
ggsave(p, filename = paste0(work.dir, "results/qcovs_boxplot.png"),
       width = 15, height = 6)
