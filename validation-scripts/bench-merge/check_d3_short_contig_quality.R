rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)
library(ggrepel)

# work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
work.dir <- "../../../revision/kamrat_ctg_align/"

min_len <- 61

evaluate_local <- function(ctg.fa, align.res) {
    perfect.align <- (align.res$qlen == align.res$qend &
                          align.res$qstart == 1 &
                          align.res$qlen == align.res$length &
                          align.res$qlen == align.res$nident &
                          align.res$pident == 100)
    return(round(sum(perfect.align) / length(ctg.fa) * 100, 2))
}

evaluate_global <- function(ctg.fa, align.res) {
    ident <- (align.res$qlen == align.res$qend &
                  align.res$qstart == 1 &
                  align.res$qlen == align.res$length &
                  align.res$qlen == align.res$nident &
                  align.res$sstart == 1 &
                  align.res$send == align.res$slen &
                  align.res$slen == align.res$qlen &
                  align.res$pident == 100)
    return(round(sum(ident) / length(ctg.fa) * 100, 2))
}

pdf(paste0(work.dir, "results/check3_cmp_short_ctg_quality.pdf"),
    width=9, height=3)
for (abd in c(1, 5)) {
    stats.res <- NULL
    for (dpt in c(1, 5, 10, 20, 30, 40, 50)) {
        # KaMRaT none
        cat(dpt, "KaMRaT none", "\n")
        ctg.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-", abd, "/depth_", dpt, "/ctg-seq.none.fa")
        align.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-", abd, "/depth_", dpt, "/ctg-aligned.none.tsv")
        ctg.fa <- readDNAStringSet(ctg.path)
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
        cat(length(ctg.fa), nrow(align.res), "\n")
        ctg.fa <- ctg.fa[width(ctg.fa) <= min_len]
        align.res <- align.res[align.res$qlen <= min_len, ]
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = "KaMRaT none",
                                      "nb.ctg" = length(ctg.fa),
                                      "ctg.median.len" = median(nchar(ctg.fa)),
                                      "perf.align" = evaluate_local(ctg.fa, align.res),
                                      "identical" = evaluate_global(ctg.fa, align.res)))
        # KaMRaT with intervention
        for (mode in c("pearson", "spearman", "mac")) {
            cat(dpt, "KaMRaT", paste(mode, 0.2, sep = ":"), "\n")
            ctg.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-", abd, "/depth_", dpt, "/ctg-seq.", mode, "_0.2", ".fa")
            align.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-", abd, "/depth_", dpt, "/ctg-aligned.", mode, "_0.2", ".tsv")
            ctg.fa <- readDNAStringSet(ctg.path)
            align.res <- read.table(align.path, header = TRUE, row.names = 1)
            cat(length(ctg.fa), nrow(align.res), "\n")
            ctg.fa <- ctg.fa[width(ctg.fa) <= min_len]
            align.res <- align.res[align.res$qlen <= min_len, ]
            stats.res <- rbind(stats.res,
                               data.frame("depth" = dpt,
                                          "mode" = paste0("KaMRaT ", mode, ":0.2"),
                                          "nb.ctg" = length(ctg.fa),
                                          "ctg.median.len" = median(nchar(ctg.fa)),
                                          "perf.align" = evaluate_local(ctg.fa, align.res),
                                          "identical" = evaluate_global(ctg.fa, align.res)))
        }
    }
    write.csv(stats.res, paste0(work.dir, "results/check3_cmp_short_ctg_quality_2-", abd, ".csv"),
              quote = FALSE)
    
    plt <- ggplot(data = stats.res, aes(x = depth, y = perf.align, color = mode)) +
        geom_line(linewidth = 1) +
        geom_point() +
        # geom_text_repel(aes(label = ctg.median.len)) +
        scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 40, 50)) +
        scale_color_manual(values = c("KaMRaT none" = "#e66101",
                                      "KaMRaT mac:0.2" = "#fdb863",
                                      "KaMRaT pearson:0.2" = "#b2abd2",
                                      "KaMRaT spearman:0.2" = "#5e3c99",
                                      "rnaSPAdes allkmers" = "#808080",
                                      "rnaSPAdes allreads" = "#000000")) +
        ylab("%perfect alignment") +
        theme_light() +
        theme(text = element_text(size = 15, family = "sans"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    plot(plt)
}
dev.off()
