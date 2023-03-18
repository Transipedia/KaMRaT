rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)
library(ggrepel)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

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

pdf(paste0(work.dir, "results/3_cmp_thresholds.pdf"),
    width=9, height=3)
stats.res <- NULL
for (thres in seq(0.1, 1, 0.1)) {
    # KaMRaT with intervention
    for (mode in c("pearson", "spearman", "mac")) {
        cat("KaMRaT", paste(mode, thres, sep = ":"), "\n")
        ctg.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-1/depth_10/ctg-seq.", mode, "_",
			   formatC(thres, digits = 1, format = "f"), ".fa")
        align.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-1/depth_10/ctg-aligned.", mode, "_",
			     formatC(thres, digits = 1, format = "f"), ".tsv")
        ctg.fa <- readDNAStringSet(ctg.path)
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
	ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
	align.res <- align.res[align.res$qlen > min_len, ]
        stats.res <- rbind(stats.res,
                           data.frame("threshold" = thres,
                                      "mode" = paste0("KaMRaT ", mode),
        			      "nb.ctg" = length(ctg.fa),
        			      "ctg.median.len" = median(nchar(ctg.fa)),
                                      "perf.align" = evaluate_local(ctg.fa, align.res),
        			      "identical" = evaluate_global(ctg.fa, align.res)))
    }
}
write.csv(stats.res, paste0(work.dir, "results/3_cmp_thresholds.csv"),
          quote = FALSE)

plt <- ggplot(data = stats.res, aes(x = threshold, y = perf.align, color = mode)) +
    geom_line(linewidth = 1) +
    geom_point() +
    geom_text_repel(aes(label = ctg.median.len)) +
    scale_x_reverse(breaks = seq(0, 1, 0.1)) +
    scale_color_manual(values = c("KaMRaT mac" = "#fdb863", 
        			  "KaMRaT pearson" = "#b2abd2",
        			  "KaMRaT spearman" = "#5e3c99")) + 
    ylab("%perfect alignment") +
    theme_light() +
    theme(text = element_text(size = 15, family = "sans"),
	  panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot(plt)
dev.off()
