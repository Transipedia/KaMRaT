rm(list = ls())

library(stringr)
library(magrittr)
library(corrplot)

cmdArgs <- commandArgs(trailingOnly = T)
rank.dir <- cmdArgs[1]
out.dir <- cmdArgs[2]

rank.dir <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/h_kamrat_rank_res"
out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/i_rank_cmp"

rank.res.path <- dir(rank.dir)
rank.res.path <- rank.res.path[str_detect(rank.res.path, pattern = "ranked-counts")]

rank.res <- lapply(rank.res.path, FUN = function(p) read.table(paste0(rank.dir, "/", p), header = T)[, 1 : 2])

rk.cmp <- NULL
for (r in rank.res) {
    m <- str_extract(names(r)[2], pattern = "relat_sd|ttest|snr|lrc|nbc|svm")
    if (m == "ttest" || m == "svm") {
        r[, paste0("rk.", m)] <- rank(r[, 2], ties.method = "average")
    } else if (m == "relat_sd" || m == "lrc" || m == "nbc") {
        r[, paste0("rk.", m)] <- rank(-r[, 2], ties.method = "average")
    } else if (m == "snr") {
        r[, paste0("rk.", m)] <- rank(-abs(r[, 2]), ties.method = "average")
    }
    if (is.null(rk.cmp)) {
        rk.cmp <- r[, c(1, 3)]
    } else {
        rk.cmp <- merge(rk.cmp, r[, c(1, 3)], by = "feature")
    }
}

cor.mat <- cor(rk.cmp[, 2 : ncol(rk.cmp)], use = "complete.obs")

png(filename = paste0(out.dir, "/corrplot_top500_among_methods.png"), width = 1100, height = 900)
par(mai = c(6, 0, 5, 3))
my.col <- colorRampPalette(colors = c("gold3", "white", "royalblue"))
corrplot(corr = cor.mat, order = "AOE", type = "upper", tl.pos = "d", tl.cex = 2, col = my.col(100), cl.cex = 3)
corrplot(corr = cor.mat, col = c("gold3", "royalblue"), number.cex = 3,
         add = TRUE, type = "lower", method = "number", order = "AOE", diag = FALSE, tl.pos = "n", cl.pos = "n")
dev.off()
