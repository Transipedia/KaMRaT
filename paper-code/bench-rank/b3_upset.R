rm (list = ls())

library(stringr)
library(ggplot2)
library(ggupset)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-rank/"
rank.dir <- paste0(work.dir, "kamrat_res")

dlist <- c("ft20000-es10-pout0", "ft20000-es1.5-pout0", "ft20000-es1.5-pout20", "ft200000-es1.5-pout20")
clist <- c("#5e3c99", "#b2abd2", "#fdb863", "#e66101")

for (i in 1 : 4) {
    d <- dlist[i]
    top.ft.list <- list()
    print(d)
    for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "bayes", "lr")) {
        f <- paste0(rank.dir, "/", d, "/ranked.", m, ".tsv")
        top.ft.list[[m]] <- as.character(read.table(f, header = T)[1 : 500, 1])
    }
    cmp.ft.df <- data.frame("feature" = unique(unname(unlist(top.ft.list))))
    print(nrow(cmp.ft.df))
    cmp.ft.df$mtd.list <- lapply(cmp.ft.df$feature, 
                                 FUN = function(f) {
                                     mtd.list <- NULL
                                     for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "bayes", "lr")) {
                                         if (f %in% top.ft.list[[m]]) {
                                             mtd.list <- c(mtd.list, m)
                                         }
                                     }
                                     return(mtd.list)
                                 })
    p <- ggplot(cmp.ft.df, aes(x = mtd.list)) +
        geom_bar(fill = clist[i]) +
        geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 7) +
        scale_x_upset() +
        ylim(c(0, 500)) +
        theme(text = element_text(size = 35, family = "Arial"),
              axis.title.x = element_blank(), 
              plot.margin = margin(l = 1.5, t = 0.5, unit = "cm")) +
        theme_combmatrix(combmatrix.label.text = element_text(size = 25, family = "Arial"),
                         combmatrix.panel.line.size = 1,
                         combmatrix.label.make_space = F)
    ggsave(paste0(work.dir, "/", d, "-upset.svg"), plot = p, width = 18, height = 6)
}
