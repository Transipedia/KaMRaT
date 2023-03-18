rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

evaluate_local <- function(kmer.fa, align.res) {
    perfect.align <- (align.res$qlen == align.res$qend &
                          align.res$qstart == 1 &
                          align.res$qlen == align.res$length &
                          align.res$qlen == align.res$nident &
                          align.res$pident == 100)
    return(round(sum(perfect.align) / length(kmer.fa) * 100, 2))
}
# /data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/check_kmer_error/err-illumina5_2-5
pdf(paste0(work.dir, "results/check2_kmer_errorrate.pdf"),
    width = 9, height = 5)
for (mode in c("err-free_1-1", "err-illumina5_2-1", "err-illumina5_2-5")) {
    stats.res <- NULL
    for (dpt in c(1, 5, 10, 20, 30, 40, 50)) {
        print(paste(mode, dpt))
        kmer.fa <- readDNAStringSet(paste0(work.dir, "check_kmer_error/", mode, "/depth_", dpt, ".fa"))
        align.path <- paste0(work.dir, "check_kmer_error/", mode, "/depth_", dpt, ".align.tsv")
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "error.mode" = mode,
                                      "nb.read" = nrow(align.res),
                                      "perf.align" = evaluate_local(kmer.fa, align.res)))
    }
    write.csv(stats.res, paste0(work.dir, "results/check2_kmer_errorrate.", mode, ".csv"),
              quote = FALSE)
    plt <- ggplot(stats.res, aes(x = depth, y = perf.align)) +
        geom_point() +
	geom_line(linewidth = 1) +
	scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 40, 50)) +
        ggtitle(mode) +
        theme_light() +
        theme(text = element_text(size = 15, family = "sans"))
    plot(plt)
}
dev.off()
