rm(list = ls())

library(tidyr)
library(ggplot2)
library(ggrepel)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

min_len <- 0


pdf(paste0(work.dir, "results/check4_ctglen_distribution.pdf"),
    width=9, height=5)
for (mode in c("pearson", "spearman", "mac")) {
    for (thres in c(0.1, 0.5, 1)) {
        # KaMRaT with intervention
        cat("KaMRaT", paste(mode, thres, sep = ":"), "\n")
        align.path <- paste0(work.dir, "kamrat_res_err-illumina5_2-1/depth_10/ctg-aligned.", mode, "_",
			     formatC(thres, digits = 1, format = "f"), ".tsv")
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
	align.res <- align.res[align.res$qlen > min_len, ]
	align.res$threshold <- thres
    	plt <- ggplot(data = align.res, aes(x = qlen)) +
            geom_histogram(binwidth = 1, fill = "gray") +
	    scale_y_log10(limits = c(1, 1.5e5)) +
            theme_light() +
	    ggtitle(paste0(mode, ":", thres)) +
            theme(text = element_text(size = 15, family = "sans"),
   	          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        plot(plt)
    }
}
dev.off()
