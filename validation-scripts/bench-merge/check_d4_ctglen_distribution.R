rm(list = ls())

library(tidyr)
library(ggplot2)
library(ggrepel)


# work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
work.dir <- "../../../revision/"

min_len <- 61

pdf(paste0(work.dir, "results/check4_ctglen_distribution_summarised.pdf"),
    width=9, height=5)
for (mode in c("pearson", "spearman", "mac")) {
    align.res <- NULL
    for (thres in seq(0.1, 1, 0.1)) {
        # KaMRaT with intervention
        cat("KaMRaT", paste(mode, thres, sep = ":"), "\n")
        align.path <- paste0(work.dir, "kamrat_ctg_align/kamrat_res_err-illumina5_2-1/depth_10/ctg-aligned.",
                             mode, "_", formatC(thres, digits = 1, format = "f"), ".tsv")
        align.res.x <- read.table(align.path, header = TRUE, row.names = 1)
        align.res.x$perf.align <- ifelse(align.res.x$qlen == align.res.x$qend &
                                             align.res.x$qstart == 1 &
                                             align.res.x$qlen == align.res.x$length &
                                             align.res.x$qlen == align.res.x$nident &
                                             align.res.x$pident == 100,
                                         yes = "aligned perfectly", no = "others")
        align.res.x$threshold <- as.character(thres)
        align.res <- rbind(align.res, align.res.x)
    }
    align.res$threshold <- factor(align.res$threshold, levels = as.character(seq(1, 0, -0.1)))
    plt <- ggplot(data = align.res) +
        geom_bin2d(aes(x = qlen, y = threshold),
                   binwidth = c(1, 0.1)) +
        facet_grid(perf.align ~ .) +
        scale_fill_gradient(trans = "log10", breaks = c(1, 1e2, 1e4, 1e6)) +
        scale_x_continuous(breaks = c(31, 61, 100, 150, 200, 250)) +
        ggtitle(paste0(mode)) +
        theme_light() +
        theme(text = element_text(size = 15, family = "sans"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    plot(plt)
}
dev.off()
