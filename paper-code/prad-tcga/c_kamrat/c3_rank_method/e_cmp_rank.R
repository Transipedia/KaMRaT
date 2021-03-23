rm(list = ls())

library(Biostrings)
library(stringr)
library(magrittr)
library(corrplot)

rk.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/c_kamratRes/c3_rank-methods/train"

for (i in 0 : 4) {
    all.rk <- NULL
    for (m in c("lrc", "snr", "ttest", "nbc", "svm")) {
        print(paste0(i, "-", m))
        f <- paste0(rk.dir, i, "/ranked-contig-counts.", m, ".fa")
        ctg.fa <- readDNAStringSet(f)      
        ctg.rk <- data.frame("contig" = as.character(ctg.fa))
        ctg.rk[, m] <- str_extract(names(ctg.fa), pattern = "\\_[0-9]+") %>%
            str_remove(pattern = "\\_") %>%
            as.numeric()
        if (is.null(all.rk)) {
            all.rk <- ctg.rk
        } else {
            all.rk <- merge(all.rk, ctg.rk, by = "contig")
        }
    }
    cor.mat <- cor(all.rk[, 2 : ncol(all.rk)], use = "complete.obs")

    png(filename = paste0(rk.dir, i, "/corrplot_among_methods.png"), width = 800, height = 800)
    par(family = "Arial")
    # corrplot(corr = cor.mat, order = "hclust", type = "", tl.pos = "d", tl.cex = 2, col = my.col(100),
             # cl.cex = 2)
    corrplot(corr = cor.mat, method = "color", 
             col = colorRampPalette(colors = c("gold3", "white", "royalblue"))(100), number.cex = 2, addCoef.col = "black",
             tl.col = "black", tl.cex = 2, tl.pos = "ld", tl.srt = 0, tl.offset = 1,
             cl.pos = "n",
             type = "lower", order = "hclust", diag = FALSE)
    dev.off()
}
