rm(list = ls())

library(Biostrings)
library(magrittr)
library(stringr)
library(tidyr)
library(UpSetR)
library(corrplot)
library(extrafont)

loadfonts(device = "postscript")

workdir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first"
outdir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/cmp_res/candidate-contigs"

my.col <- colorRampPalette(colors = c("white", "darkred"))

mlist <- c("ttest", "snr", "lrc", "nbc")
for (d in c("risk", "relapse-new")) {
    # for (i in 0 : 4) {
    i = 2
        print(paste0(d, "-", i))
        mtd.topctg <- lapply(mlist, FUN = function(m) {
            p <- paste0(workdir, "/", d, "/train", i, "/", m, "/signatures-ridge.fa")
            return(readDNAStringSet(p) %>% as.character())})
        names(mtd.topctg) <- mlist
        
        png(paste0(outdir, "/", d, "-upset-candidates-train", i, ".png"), width = 900, height = 500, family = "Arial")
        upset(data = fromList(mtd.topctg), order.by = "freq", nintersects = 15,
              text.scale = 3, set_size.show = T) %>% show()
        dev.off()
        
        mtd.rkctg <- lapply(mlist, FUN = function(m) {
            p <- paste0(workdir, "/", d, "/train", i, "/", m, "/ranked-contigs.fa")
            fa <- readDNAStringSet(p)
            df <- data.frame("ctg" = as.character(fa), 
                             "rk.num" = names(fa) %>% str_extract(pattern = "[0-9]+") %>% as.numeric(),
                             "rk.mtd" = m)}) %>% do.call(what = rbind)
        mtd.rkctg.wider <- pivot_wider(mtd.rkctg, id_cols = ctg, names_from = "rk.mtd", values_from = "rk.num")
        cor.mat <- cor(mtd.rkctg.wider[, -1], use = "complete.obs")
        png(paste0(outdir, "/", d, "-corr-all-train", i, ".png"), width = 900, height = 500, family = "Arial")
        corrplot(corr = cor.mat, order = "hclust", type = "upper", method = "number", col = my.col(100), cl.pos = "n",
                 tl.pos = "d", tl.cex = 3, tl.col = "black", number.cex = 3) %>% show()
        dev.off()
        
        hc <- hclust(as.dist(1 - cor.mat))
        png(paste0(outdir, "/", d, "-hclust-all-train", i, ".png"), width = 900, height = 500, family = "Arial")
        plot(hc, xlab = "Height", lwd = 12, hang = -1, cex = 3)
        dev.off()
    # }
}
