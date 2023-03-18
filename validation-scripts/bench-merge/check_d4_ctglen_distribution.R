rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)
library(patchwork)

work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
#work.dir <- "../../../revision/"


pdf(paste0(work.dir, "results/check4_ctglen_distribution_summarised.pdf"),
    width=15, height=10)
align.res <- NULL
for (mode in c("pearson", "spearman", "mac")) {
    for (thres in seq(0.1, 0.9, 0.2)) {
        # KaMRaT with intervention
        cat("KaMRaT", paste(mode, thres, sep = ":"), "\n")
        align.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-1/depth_10/ctg-aligned.",
                             mode, "_", formatC(thres, digits = 1, format = "f"), ".tsv")
        align.res.x <- read.table(align.path, header = TRUE, row.names = 1)
        align.res.x$perf.align <- ifelse(align.res.x$qlen == align.res.x$qend &
                                             align.res.x$qstart == 1 &
                                             align.res.x$qlen == align.res.x$length &
                                             align.res.x$qlen == align.res.x$nident &
                                             align.res.x$pident == 100,
                                         yes = "yes", no = "no")
        align.res.x <- align.res.x[, c("qlen", "perf.align")]
        fa.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-1/depth_10/ctg-seq.",
                          mode, "_", formatC(thres, digits = 1, format = "f"), ".fa")
        ctg.fa <- readDNAStringSet(fa.path)
        ctg.fa <- ctg.fa[!(names(ctg.fa) %in% rownames(align.res.x))]
        unaligned.res <- data.frame("qlen" = nchar(ctg.fa),
                                    "perf.align" = "no")
        rownames(unaligned.res) <- names(ctg.fa)
        align.res.x <- rbind(align.res.x, unaligned.res)
        align.res.x$threshold <- as.character(thres)
        align.res.x$mode <- mode
        align.res <- rbind(align.res, align.res.x)
    }
}
align.res$threshold <- factor(align.res$threshold, levels = as.character(seq(0.9, 0.1, -0.2)))
align.res$panel <- factor(paste(align.res$mode, align.res$perf.align, sep = "-"),
                          levels = c("pearson-yes", "spearman-yes", "mac-yes",
                                     "pearson-no", "spearman-no", "mac-no"))
plt <- ggplot(data = align.res) +
    geom_bin2d(aes(x = qlen, y = threshold), binwidth = c(1, 0.1)) +
    facet_grid(panel ~ .) +
    scale_fill_gradientn(trans = "log10", breaks = c(1, 1e2, 1e4, 1e6),
                         colours = c("gray95", "gray50", "black")) +
    scale_x_continuous(breaks = c(31, 61, 100, 150, 200, 250), limits = c(0, 290)) +
    theme_bw() +
    theme(text = element_text(size = 15, family = "sans"),
          axis.title.x = element_blank(), strip.text = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot(plt)
dev.off()
