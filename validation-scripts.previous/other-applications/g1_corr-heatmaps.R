rm (list = ls())

library(magrittr)
library(tidyr)
library(pheatmap)
library(ggwordcloud)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"

smp.info <- read.table(paste0(work.dir, "matrices/all-samples.tsv"), header = F, row.names = 1)
names(smp.info) <- "condition"

nb.cell <- read.table("/home/haoliang.xue/media/user/development/KaMRaT/paper-code/real-data/sample-conditions/Macrophages.M2.tsv",
                      header = F, row.names = 1)
names(nb.cell) <- "MP.M2.abundance"

smp.annot <- cbind(smp.info, nb.cell)

smp.annot <- smp.annot[order(smp.annot$MP.M2.abundance), ]

print(all(nb.cell[rownames(smp.annot), "MP.M2.abundance"] == smp.annot$MP.M2.abundance))
print(all(nb.cell[rownames(smp.annot), "condition"] == smp.info$condition))

for (f in c("pearson", "spearman")) {
    print(f)
    kamrat.res <- read.table(paste0(work.dir, "kamrat_res/all-samples/corr-rank-merge/top-contig-counts-", f, ".tsv"), 
                             header = T)
    annot.res <- read.table(paste0(work.dir, "kamrat_res/all-samples/corr-rank-merge/merged-annotation-", f, ".tsv"),
                            header = T)
    res <- merge(kamrat.res, annot.res)
    res$annot <- paste0(as.character(res$gene_symbol), "-", "ctg", 1 : nrow(res))
    
    top.res <- res[order(abs(res$rep.value), decreasing = T)[1 : 10], 
                   c("annot", "rep.value", rownames(smp.annot))]
    top.res <- pivot_longer(top.res, cols = c(-annot, -rep.value), names_to = "sample", values_to = "count")
    top.res <- merge(top.res, nb.cell, by.x = "sample", by.y = "row.names")
    
    p <- ggplot(top.res) +
        geom_point(aes(x = count, y = MP.M2.abundance)) +
        facet_wrap(annot ~ ., nrow = 2, scales = "free") +
        theme(text = element_text(size = 15, family = "Arial"))
    ggsave(p, filename = paste0(work.dir, "kamrat_res/all-samples/corr-rank-merge/corr-scatter-", f, ".svg"),
           width = 16, height = 9)
    
    rownames(res) <- res$annot
    res <- res[res$nb.merged.kmer > 30, ]
    print(as.character(unique(res$gene_symbol)))
    res <- res[, rownames(smp.annot)]
    pheatmap(res,
             cluster_cols = F,
             cluster_rows = T,
             color = colorRampPalette(c("white", "gray", "black"), bias = 3)(100),
             border_color = NA,
             annotation_col = smp.annot[c("MP.M2.abundance", "condition")],
             annotation_colors = list(condition = c("normal" = "royalblue", "tumor" = "gold3"),
                                      MP.M2.abundance = colorRampPalette(c("white", "brown4"))(100)),
             show_colnames = FALSE,
             show_rownames = TRUE,
             fontsize = 10,
             cellheight = 10,
             cellwidth = 3,
             width = 20,
             filename = paste0(work.dir, "kamrat_res/all-samples/corr-rank-merge/corr-heatmap-", f, ".tiff"))
    
}
