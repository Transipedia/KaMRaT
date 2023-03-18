rm(list = ls())

library(Biostrings)
library(stringr)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-merge/"

merge.summary <- NULL
for (m in c("none", "mac", "pearson", "spearman")) {
    for (p in seq(40, 100, by = 10)) {
        print(paste0(m, "-", p))
        merged.ctg.nb <- readDNAStringSet(paste0(work.dir, "blastn_res/", p, "pct/ctg-seq.", m, 
                                                 ifelse(m == "none", no = "_0.2.fa", yes = ".fa"))) %>% length()
        merged.ctg.align <- read.table(paste0(work.dir, "blastn_res/", p, "pct/ctg-aligned.", m, 
                                              ifelse(m == "none", no = "_0.2.tsv", yes = ".tsv")), header = T)
        perf.ctg.nb <- merged.ctg.align[merged.ctg.align$qlen == merged.ctg.align$qend &
                                            merged.ctg.align$qstart == 1 &
                                            merged.ctg.align$qlen == merged.ctg.align$length &
                                            merged.ctg.align$qlen == merged.ctg.align$nident &
                                            merged.ctg.align$pident == 100, ] %>%
            nrow()
        merge.summary <- rbind(merge.summary,
                               data.frame("rk.mtd" = str_extract(m, pattern = "none|mac|pearson|spearman"),
                                          "miss.perc" = 100 - p,
                                          "perf.ratio" = perf.ctg.nb / merged.ctg.nb))        
    }
}

p <- ggplot(merge.summary, aes(x = miss.perc, y = perf.ratio * 100, color = rk.mtd)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("none" = "#e66101", "mac" = "#fdb863", "pearson" = "#b2abd2", "spearman" = "#5e3c99")) +
    scale_x_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%")) +
    scale_y_continuous(breaks = seq(100, 78, by = -2), labels = paste0(seq(100, 78, by = -2), "%")) +
    xlab("percent of missing k-mers for extension") +
    ylab("perfect alignment ratio") +
    theme(text = element_text(size = 30, family = "Arial"), legend.position = "none")
ggsave(paste0(work.dir, "/merge-perfectAlignRatio.svg"), plot = p, width = 12, height = 6)
