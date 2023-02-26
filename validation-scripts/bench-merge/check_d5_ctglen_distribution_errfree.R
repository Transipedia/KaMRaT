rm(list = ls())

library(tidyr)
library(ggplot2)
library(ggrepel)


# work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
work.dir <- "../../../revision/"

min_len <- 61

pdf(paste0(work.dir, "results/check5_ctglen_distribution_errfree.pdf"),
    width=9, height=5)
for (method in c("pearson", "spearman", "mac")) {
    align.res <- NULL
    # KaMRaT with intervention
    align.path <- paste0(work.dir, "kamrat_ctg_align/kamrat_res_err-free_1-1/depth_10/ctg-aligned.",
                         method, "_0.2.tsv")
    align.res <- read.table(align.path, header = TRUE, row.names = 1)
    align.res$perf.align <- ifelse(align.res$qlen == align.res$qend &
                                       align.res$qstart == 1 &
                                       align.res$qlen == align.res$length &
                                       align.res$qlen == align.res$nident &
                                       align.res$pident == 100,
                                     yes = "aligned perfectly", no = "others")
    plt <- ggplot(data = align.res, aes(x = qlen)) +
        geom_histogram(binwidth = 1) +
        facet_grid(perf.align ~ ., scales = "free_y") +
        scale_x_continuous(limits = c(31, 300),
                           breaks = c(31, 61, 100, 150, 200, 250, 300)) +
        ggtitle(paste0(method, ":0.2")) +
        theme_light() +
        theme(text = element_text(size = 15, family = "sans"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    plot(plt)
}
dev.off()
