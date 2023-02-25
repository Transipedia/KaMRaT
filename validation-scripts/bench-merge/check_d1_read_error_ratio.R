rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

evaluate_local <- function(align.res) {
    perfect.align <- (align.res$qlen == align.res$qend &
                          align.res$qstart == 1 &
                          align.res$qlen == align.res$length &
                          align.res$qlen == align.res$nident &
                          align.res$pident == 100)
    return(round(sum(perfect.align) / nrow(align.res) * 100, 2))
}

pdf(paste0(work.dir, "results/check1_read_errorrate.pdf"),
    width = 9, height = 5)
for (mode in c("err-free", "err-illumina5")) {
    stats.res <- NULL
    for (dpt in as.character(c(1, 5, 10, 20, 30, 40, 50))) {
        print(paste(mode, dpt))
        depth.dir <- paste0(work.dir, "check_read_error/", mode, "/depth_", dpt)
        for (f in dir(depth.dir)) {
            align.path <- paste0(depth.dir, "/", f)
            align.res <- read.table(align.path, header = TRUE, row.names = 1)
            stats.res <- rbind(stats.res,
                               data.frame("depth" = dpt,
                                          "read.file" = strsplit(f, "\\.")[[1]][1],
                                          "error.mode" = mode,
                                          "nb.read" = nrow(align.res),
                                          "perf.align" = evaluate_local(align.res)))
        }
    }
    write.csv(stats.res, paste0(work.dir, "results/check1_read_errorrate.", mode, ".csv"),
              quote = FALSE)
    plt <- ggplot(stats.res, aes(x = depth, y = perf.align)) +
        geom_boxplot(fill = "gray") +
        ggtitle(mode) +
        geom_jitter(height = 0) +
        theme_light() +
        theme(text = element_text(size = 15, family = "sans"))
    plot(plt)
}
dev.off()
