rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)
library(patchwork)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge-above5500/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

min_len <- 61

evaluate_perfalign <- function(ctg.fa, align.res) {
    perfect.align <- (align.res$qlen == align.res$qend &
                      align.res$qstart == 1 &
                      align.res$qlen == align.res$length &
                      align.res$qlen == align.res$nident &
                      align.res$pident == 100 & 
		      align.res$qlen > min_len)
    return(sum(perfect.align) / sum(width(ctg.fa) > min_len) * 100)
}

evaluate_compalign <- function(ctg.fa, align.res) {
    return(sum(align.res$qcovs == 100 & align.res$qlen > min_len) / sum(width(ctg.fa) > min_len) * 100)
}

abd <- 1
stats.res <- NULL
for (dpt in c(0.05, 0.2, 0.4, 0.6, 0.8, 1)) {
    # KaMRaT none
    cat(dpt, "KaMRaT none", ":")
    ctg.path <- paste0(work.dir, "kamrat_res_err-free_1-", abd, "/depth_", dpt, "/ctg-seq.none.fa")
    align.path <- paste0(work.dir, "kamrat_res_err-free_1-", abd, "/depth_", dpt, "/ctg-aligned.none.tsv")
    ctg.fa <- readDNAStringSet(ctg.path)
    align.res <- read.table(align.path, header = TRUE, row.names = 1)
    ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
    align.res <- align.res[align.res$qlen > min_len, ]
    cat("\t", length(ctg.fa), nrow(align.res), "\n")
    stats.res <- rbind(stats.res,
                       data.frame("depth" = dpt,
                                  "mode" = "KaMRaT none",
				  "nb.ctg" = length(ctg.fa),
				  "ctg.median.len" = median(nchar(ctg.fa)),
				  "perf.align" = evaluate_perfalign(ctg.fa, align.res),
				  "comp.align" = evaluate_compalign(ctg.fa, align.res)))
    # KaMRaT with intervention
    for (mode in c("pearson", "spearman", "mac")) {
        cat(dpt, "KaMRaT", paste(mode, 0.2, sep = ":"), ":")
        ctg.path <- paste0(work.dir, "kamrat_res_err-free_1-", abd, "/depth_", dpt, "/ctg-seq.", mode, "_0.2", ".fa")
        align.path <- paste0(work.dir, "kamrat_res_err-free_1-", abd, "/depth_", dpt, "/ctg-aligned.", mode, "_0.2", ".tsv")
        ctg.fa <- readDNAStringSet(ctg.path)
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
        ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
        align.res <- align.res[align.res$qlen > min_len, ]
        cat("\t", length(ctg.fa), nrow(align.res), "\n")
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = paste0("KaMRaT ", mode, ":0.2"),
				      "nb.ctg" = length(ctg.fa),
				      "ctg.median.len" = median(nchar(ctg.fa)),
				      "perf.align" = evaluate_perfalign(ctg.fa, align.res),
				      "comp.align" = evaluate_compalign(ctg.fa, align.res)))
    }
    # rnaSPAdes
    for (mode in c("allreads", paste0("allkmers-1-", abd))) {
        cat(dpt, "SPAdes", mode, ":")
        ctg.path <- paste0(work.dir, "spades_res/err-free/depth_", dpt, "/", mode, "/transcripts.fasta")
        align.path <- paste0(work.dir, "spades_res/err-free/depth_", dpt, "/", mode, "/blastn_align.tsv")
        ctg.fa <- readDNAStringSet(ctg.path)
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
        cat("\t", length(ctg.fa), nrow(align.res), "\n")
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = paste0("rnaSPAdes ", strsplit(mode, split = "-")[[1]][1]),
				      "nb.ctg" = length(ctg.fa),
				      "ctg.median.len" = median(nchar(ctg.fa)),
				      "perf.align" = evaluate_perfalign(ctg.fa, align.res),
				      "comp.align" = evaluate_compalign(ctg.fa, align.res)))
    }
}
write.csv(stats.res, paste0(work.dir, "results/2_newcmp_with_without_intervention_1-", abd, ".csv"),
          quote = FALSE)

pdf(paste0(work.dir, "results/2_newcmp_with_without_intervention.pdf"),
    width=9, height=7)
plt1 <- ggplot(data = stats.res, aes(x = depth, y = perf.align, color = mode)) +
    geom_line(linewidth = 1) +
    geom_point() +
    scale_x_continuous(breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "KaMRaT mac:0.2" = "#fdb863",
				  "KaMRaT pearson:0.2" = "#b2abd2",
				  "KaMRaT spearman:0.2" = "#5e3c99",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylim(c(70, 100)) +
    ylab("%perfect alignment") +
    theme_light() +
    theme(text = element_text(size = 15, family = "sans"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 <- ggplot(data = stats.res, aes(x = depth, y = comp.align, color = mode)) +
    geom_line(linewidth = 1) +
    geom_point() +
    scale_x_continuous(breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "KaMRaT mac:0.2" = "#fdb863",
				  "KaMRaT pearson:0.2" = "#b2abd2",
				  "KaMRaT spearman:0.2" = "#5e3c99",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylab("%complete alignment") +
    theme_light() +
    theme(text = element_text(size = 15, family = "sans"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot(plt1/plt2)
dev.off()
