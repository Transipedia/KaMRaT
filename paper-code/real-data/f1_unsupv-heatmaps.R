rm (list = ls())

library(magrittr)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(pheatmap)
# library(RColorBrewer)

# work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"
work.dir <- "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"

smp.info <- read.table(paste0(work.dir, "matrices/all-samples.tsv"), header = F, row.names = 1)
names(smp.info) <- "condition"

for (f in c("entropy", "sd", "rsd1", "rsd2", "rsd3")) {
    print(f)
    kamrat.res <- read.table(paste0(work.dir, "kamrat_res/all-samples/unsupv-rank-merge/contig-counts-", f, ".tsv"), 
                             header = T, row.names = 1)[, c(-1, -2, -3)] %>% t()
    print(ncol(kamrat.res))
    kamrat.dend <- dist(kamrat.res) %>% hclust() %>% as.dendrogram()
    p <- kamrat.dend %>%
        set("branches_lwd", 0.5) %>%
        set("leaves_col",
            value = ifelse(smp.info[labels(kamrat.dend), "condition"] == "tumor",
                           yes = "gold3", no = "royalblue")) %>%
        set("leaves_pch", 16) %>%
        set("leaves_cex", 3) %>%
        as.ggdend() %>%
        ggplot(labels = FALSE) + scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
    ggsave(p, filename = paste0(work.dir, "kamrat_res/all-samples/unsupv-rank-merge/unsupv-dendrogram-", f, ".svg"),
           width = 9, height = 9)

    pheatmap(t(kamrat.res),
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             color = colorRampPalette(c("white", "black"), bias = 5)(50),
             border_color = NA,
             annotation_col = smp.info[1],
             annotation_colors = list(condition = c("normal" = "royalblue", "tumor" = "gold3")),
             show_colnames = FALSE,
             show_rownames = FALSE,
             fontsize = 18,
             cellwidth = 3,
             width = 20,
             filename = paste0(work.dir, "kamrat_res/all-samples/unsupv-rank-merge/unsupv-heatmap-", f, ".tiff"))
}
