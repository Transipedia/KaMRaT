rm(list = ls())

library(Biostrings)
library(stringr)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-merge/"

merge.summary <- NULL
for (m in c("mac", "pearson", "spearman")) {
    for (thres in seq(0.1, 1, by = 0.1)) {
        print(paste0(m, "-", sprintf("%.1f", thres)))
        merged.ctg <- readDNAStringSet(paste0(work.dir, "blastn_res/40pct/ctg-seq.", m, "_", 
                                              sprintf("%.1f", thres), ".fa"))
        merged.ctg.align <- read.table(paste0(work.dir, "blastn_res/40pct/ctg-aligned.", m, "_", 
                                              sprintf("%.1f", thres), ".tsv"), header = T)
        perf.ctg.nb <- merged.ctg.align[merged.ctg.align$qlen == merged.ctg.align$qend &
                                            merged.ctg.align$qstart == 1 &
                                            merged.ctg.align$qlen == merged.ctg.align$length &
                                            merged.ctg.align$qlen == merged.ctg.align$nident &
                                            merged.ctg.align$pident == 100, ] %>%
            nrow()
        merge.summary <- rbind(merge.summary,
                               data.frame("rk.mtd" = m,
                                          "threshold" = sprintf("%.1f", thres),
                                          "mean.len" = mean(nchar(merged.ctg)),
                                          "sd.len" = sd(nchar(merged.ctg)),
                                          "median.len" = median(nchar(merged.ctg)),
                                          "perf.ratio" = perf.ctg.nb / length(merged.ctg)))        
    }
}

p <- ggplot() +
    geom_line(data = merge.summary, aes(x = median.len, y = perf.ratio * 100, color = rk.mtd), size = 1) +
    geom_point(data = merge.summary, aes(x = median.len, y = perf.ratio * 100, color = rk.mtd, 
                                         size = as.character(threshold))) +
    scale_color_manual(values = c("none" = "#e66101", "mac" = "#fdb863", 
                                  "pearson" = "#b2abd2", "spearman" = "#5e3c99")) +
    scale_x_continuous(breaks = seq(50, 100, by = 5)) +
    scale_y_continuous(breaks = seq(78, 100, by = 2),
                       labels = paste0(seq(78, 100, by = 2), "%")) +
    xlab("contig median length (nt)") +
    ylab("perfect alignment ratio") +
    theme(text = element_text(size = 30, family = "Arial"), legend.position = "none")
ggsave(paste0(work.dir, "/merge-threshold.svg"), plot = p, width = 12, height = 6)
