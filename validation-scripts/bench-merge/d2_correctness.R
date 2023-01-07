rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)

work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

evaluate_align <- function(ctg.path, align.path) {
    ctg.fa <- readDNAStringSet(ctg.path)
    align.res <- read.table(align.path, header = TRUE, row.names = 1)
    if (!all(nchar(ctg.fa) == align.res$qlen[rownames(ctg.fa)])) {
        stop(paste(ctg.path, align.path, sep = "\n"))
    }
    align.res$perfect.align <- (align.res$qlen == align.res$qend &
                                align.res$qstart == 1 &
                                align.res$qlen == align.res$length &
                                align.res$qlen == align.res$nident &
                                align.res$pident == 100)
    return(sum(align.res$perfect.align) / length(ctg.fa) * 100)
}

stats.res <- NULL
for (dpt in c(0.1, 0.2, 0.5, 1, 2, 5, 10)) {
    # KaMRaT none
    cat(dpt, "KaMRaT none", "\n")
    ctg.path <- paste0(work.dir, "kamrat_res/depth_", dpt, "/ctg-seq.none.fa")
    align.path <- paste0(work.dir, "kamrat_res/depth_", dpt, "/ctg-aligned.none.tsv")
    stats.res <- rbind(stats.res,
                       data.frame("depth" = dpt,
                                  "mode" = "KaMRaT none",
                                  "perf.align" = evaluate_align(ctg.path, align.path)))
    # KaMRaT with intervention
    for (mode in c("pearson", "spearman", "mac")) {
        cat(dpt, "KaMRaT", paste(mode, 0.2, sep = ":"), "\n")
        ctg.path <- paste0(work.dir, "kamrat_res/depth_", dpt, "/ctg-seq.", mode, "_0.2", ".fa")
        align.path <- paste0(work.dir, "kamrat_res/depth_", dpt, "/ctg-aligned.", mode, "_0.2", ".tsv")
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = paste0("KaMRaT ", mode, ":0.2"),
                                      "perf.align" = evaluate_align(ctg.path, align.path)))
    }
    # rnaSPAdes
    for (mode in c("allreads", "allkmers")) {
        cat(dpt, "SPAdes", mode, "\n")
        ctg.path <- paste0(work.dir, "spades_res/depth_", dpt, "/", mode, "/transcripts.fasta")
        align.path <- paste0(work.dir, "spades_res/depth_", dpt, "/", mode, "/blastn_align.tsv")
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = paste("rnaSPAdes", mode),
                                      "perf.align" = evaluate_align(ctg.path, align.path)))
    }
}

write.csv(stats.res, paste0(work.dir, "results/correctness.csv"), quote = FALSE)

pdf(paste0(work.dir, "results/correctness.pdf"), width=16, height=10)
ggplot(data = stats.res) +
    geom_line(aes(x = depth, y = 100 - perf.align, color = mode)) +
    geom_point(aes(x = depth, y = 100 - perf.align, color = mode)) +
    scale_x_log10(breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "KaMRaT mac:0.2" = "#fdb863", 
				  "KaMRaT pearson:0.2" = "#b2abd2",
				  "KaMRaT spearman:0.2" = "#5e3c99",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylab("mis-extension/-assenmbly %") +
    theme_light() +
    theme(text = element_text(size = 30, family = "sans"))
ggplot(data = stats.res) +
    geom_line(aes(x = depth, y = 100 - perf.align, color = mode)) +
    geom_point(aes(x = depth, y = 100 - perf.align, color = mode)) +
    scale_x_log10(breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    scale_y_log10() +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "KaMRaT mac:0.2" = "#fdb863", 
				  "KaMRaT pearson:0.2" = "#b2abd2",
				  "KaMRaT spearman:0.2" = "#5e3c99",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylab("mis-extension/-assenmbly %") +
    theme_light() +
    theme(text = element_text(size = 30, family = "sans"))
dev.off()
