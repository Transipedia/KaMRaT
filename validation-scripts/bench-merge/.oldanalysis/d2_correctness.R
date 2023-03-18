rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)

work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

evaluate_align <- function(ctg.path, align.path) {
    ctg.fa <- readDNAStringSet(ctg.path)
    # ctg.fa <- ctg.fa[width(ctg.fa) > 61]
    align.res <- read.table(align.path, header = TRUE, row.names = 1)
    #align.res <- align.res[align.res$qlen > 61, ]
    if (!all(nchar(ctg.fa) == align.res$qlen[rownames(ctg.fa)])) {
        stop(paste(ctg.path, align.path, sep = "\n"))
    }
    perfect.align <- (align.res$qlen == align.res$qend &
                      align.res$qstart == 1 &
                      align.res$qlen == align.res$length &
                      align.res$qlen == align.res$nident &
                      align.res$pident == 100)
    return(sum(perfect.align) / length(ctg.fa) * 100)
}

stats.res <- NULL
for (toppct in seq(0.1, 1, 0.1)) {
    # KaMRaT none
    cat(toppct, "KaMRaT none", "\n")
    suffix <- ifelse(toppct == 1, yes = "", no = paste0(".topkmers-ttest.padj-", toppct))
    ctg.path <- paste0(work.dir, "kamrat_res/depth_10/ctg-seq.none", suffix, ".fa")
    align.path <- paste0(work.dir, "kamrat_res/depth_10/ctg-aligned.none", suffix, ".tsv")
    stats.res <- rbind(stats.res,
                       data.frame("top.pct" = toppct,
                                  "mode" = "KaMRaT none",
                                  "perf.align" = evaluate_align(ctg.path, align.path)))
    # KaMRaT with intervention
    for (mode in c("pearson", "spearman", "mac")) {
        cat(toppct, "KaMRaT", paste(mode, 0.2, sep = ":"), "\n")
        ctg.path <- paste0(work.dir, "kamrat_res/depth_10/ctg-seq.", mode, "_0.2", suffix, ".fa")
        align.path <- paste0(work.dir, "kamrat_res/depth_10/ctg-aligned.", mode, "_0.2", suffix, ".tsv")
        stats.res <- rbind(stats.res,
                           data.frame("top.pct" = toppct,
                                      "mode" = paste0("KaMRaT ", mode, ":0.2"),
                                      "perf.align" = evaluate_align(ctg.path, align.path)))
    }
}

write.csv(stats.res, paste0(work.dir, "results/kamrat_correctness_ttest.padj.csv"), quote = FALSE)

pdf(paste0(work.dir, "results/kamrat_correctness_ttest.padj.pdf"), width=16, height=10)
ggplot(data = stats.res) +
    geom_line(aes(x = top.pct, y = perf.align, color = mode)) +
    geom_point(aes(x = top.pct, y = perf.align, color = mode)) +
    scale_x_reverse(breaks = seq(0, 1, 0.1)) +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "KaMRaT mac:0.2" = "#fdb863", 
				  "KaMRaT pearson:0.2" = "#b2abd2",
				  "KaMRaT spearman:0.2" = "#5e3c99",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylab("%contig perfectly aligned") +
    theme_light() +
    theme(text = element_text(size = 30, family = "sans"))
dev.off()
