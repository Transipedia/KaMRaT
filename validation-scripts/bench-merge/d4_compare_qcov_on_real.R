rm(list = ls())

setwd("~/Projects/working-bench/KaMRaT/validation-scripts/bench-merge")
library(ggpubr)

read_align <- function(merge_mode) {
    align.res <- read.table(paste0("../../../revision/merge-on-real/PRADtcga/contig-alignment.",
                                   merge_mode, ".tsv"),
                            sep="\t", header = TRUE, stringsAsFactors=F)
    align.res$source <- merge_mode
    return(align.res)
}
align <- list()
for (m in c("none", "pearson", "spearman", "mac")) {
    print(m)
    align[[m]] <- read_align(m)
}

pval.df <- NULL
for(merge_mode1 in c("none", "pearson", "spearman", "mac")) {
    for (merge_mode2 in c("none", "pearson", "spearman", "mac")) {
        pval <- wilcox.test(unlist(align[[merge_mode1]]["qcovs"]),
                            unlist(align[[merge_mode2]]["qcovs"]),
                            alternative = "two.sided")$p.value
        pval.df <- rbind(pval.df,
                         data.frame("merge.mode.1" = merge_mode1,
                                    "merge.mode.2" = merge_mode2,
                                    "pvalue.nominal" = pval))
    }
}
pval.df$pvalue.adj <- p.adjust(pval.df$pvalue.nominal, method = "bonferroni")
pval.df$significant <- pval.df$pvalue.adj < 0.05
pval.df <- pval.df[order(pval.df$pvalue.adj, pval.df$pvalue.nominal, decreasing = FALSE), ]
write.csv(pval.df, file = "../../../revision/merge-on-real/PRADtcga/qcovs_wx-test.csv",
          quote = FALSE, row.names = FALSE)

p <- ggboxplot(rbind(align[["none"]], align[["mac"]], align[["pearson"]], align[["spearman"]]),
               x = "source", y = "qcovs", color = "source",
               add = "jitter", size = 1, notch = TRUE) + 
    scale_color_manual(values = c("none" = "#e66101",
                                  "mac" = "#fdb863",
                                  "pearson" = "#b2abd2",
                                  "spearman" = "#5e3c99"))
ggsave(p, filename = "../../../revision/merge-on-real/PRADtcga/qcovs_boxplot.png",
       width = 15, height = 6)
