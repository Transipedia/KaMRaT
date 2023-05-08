rm(list = ls())

library(Biostrings)
library(magrittr)
library(UpSetR)
library(ggplot2)
library(patchwork)

ctg.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/kamrat_on_varmat/apply-on-varmat/140-smp/"

read_align <- function(merge_mode) {
    align.res <- read.table(paste0(ctg.dir, "contig-alignment.", merge_mode, ".tsv"),
                            sep="\t", header = TRUE, stringsAsFactors=F, row.names = 1)
    align.fa <- readDNAStringSet(paste0(ctg.dir, "contig-seq.", merge_mode, ".fa"))
    align.rc.fa <- reverseComplement(align.fa)
    align.res$seq1 <- as.character(align.fa[rownames(align.res)])
    align.res$seq2 <- as.character(align.rc.fa[rownames(align.res)])
    align.res$seq <- apply(align.res, MARGIN = 1,
                           FUN = function(line) min(line["seq1"], line["seq2"])) %>%
        unlist()
    align.res$source <- merge_mode
    return(align.res)
}
align <- list()
for (m in c("none", "pearson", "spearman", "mac")) {
    print(m)
    align[[m]] <- read_align(m)
    align[[m]]$my_pident <- (align[[m]]$nident / align[[m]]$qlen * 100)
    cat(m, median(align[[m]]$qlen), nrow(align[[m]]), "\n")
    align[[m]] <- align[[m]][align[[m]]$qlen > 61, ]
}

ctg.shared <- intersect(x = align[["none"]]$seq, y = align[["mac"]]$seq)
ctg.shared <- intersect(x = ctg.shared, y = align[["pearson"]]$seq)
ctg.shared <- intersect(x = ctg.shared, y = align[["spearman"]]$seq)

for (m in c("none", "pearson", "spearman", "mac")) {
    print(m)
    align[[m]] <- align[[m]][!(align[[m]]$seq %in% ctg.shared), ]
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
pval.df$significant <- (pval.df$pvalue.adj < 0.05)
pval.df <- pval.df[order(pval.df$pvalue.adj, pval.df$pvalue.nominal,
                         decreasing = FALSE), ]

all.df <- rbind(align[["none"]], align[["mac"]], align[["pearson"]], align[["spearman"]])
all.df$source <- factor(all.df$source, levels = c("none", "mac", "pearson", "spearman"))
p1 <- ggplot(data = all.df, aes(x = source, y = qcovs, fill = source)) +
    geom_boxplot(outlier.size = 1, outlier.alpha = 0.1) +
    scale_fill_manual(values = c("none" = "#e66101",
                                 "mac" = "#fdb863",
                                 "pearson" = "#b2abd2",
                                 "spearman" = "#5e3c99")) +
    theme_bw() +
    theme(text = element_blank(), legend.position = "none",
          panel.grid = element_blank())
p2 <- ggplot(data = all.df, aes(x = source, y = my_pident, fill = source)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("none" = "#e66101",
                                 "mac" = "#fdb863",
                                 "pearson" = "#b2abd2",
                                 "spearman" = "#5e3c99")) +
    theme_bw() +
    theme(text = element_text(size = 15), legend.position = "none")
ggsave(p1, filename = paste0(ctg.dir, "qcovs_pidents_boxplots_notshared.png"),
       width = 9, height = 3.5, dpi = 600)
