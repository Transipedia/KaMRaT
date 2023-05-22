rm (list = ls())

library(tidyr)
library(magrittr)
library(stringr)
library(foreach)
library(doParallel)

registerDoParallel(makeCluster(6))

work.dir <- "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/"
# work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/"
out.dir <- paste0(work.dir, "RealData/")

for (d in c("PRADtcga", "LUADseo")) {
    print(d)
    
    ft.before <- list()
    for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes")) {
        print(m)
        ft.before[[m]] <- read.table(paste0(work.dir, d, "/kamrat_res/CV0/merge-rank/", 
                                            m, "/top-merged-contig-counts.tsv"), header = T)$contig %>%
            as.character()
    }
    cmp.ft.df.before <- foreach(f = unique(unname(unlist(ft.before))), .combine = rbind) %dopar% {
        library(stringr)
        
        mtd.list <- NULL
        for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes")) {
            if (f %in% ft.before[[m]]) {
                mtd.list <- c(mtd.list, m)
            }
        }
        data.frame("feature" = f, "mtd" = mtd.list)
    }
    cmp.ft.df.before <- aggregate(cmp.ft.df.before$mtd, by = list(cmp.ft.df.before$feature), FUN = paste0)
    names(cmp.ft.df.before) <- c("feature", "mtd.list")
    
    ft.after <- list()
    for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes")) {
        print(m)
        ft.after[[m]] <- read.table(paste0(work.dir, "/", d, "/kamrat_res/CV0/merge-rank/", 
                                           m, "/randomforest_models/0_RF.pickle.features"), header = F)$V1 %>%
            as.character()
    }
    cmp.ft.df.after <- foreach(f = unique(unname(unlist(ft.after))), .combine = rbind) %dopar% {
        library(stringr)
        
        mtd.list <- NULL
        for (m in c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes")) {
            if (f %in% ft.after[[m]]) {
                mtd.list <- c(mtd.list, m)
            }
        }
        data.frame("feature" = f, "mtd" = mtd.list)
    }
    cmp.ft.df.after <- aggregate(cmp.ft.df.after$mtd, by = list(cmp.ft.df.after$feature), FUN = paste0)
    names(cmp.ft.df.after) <- c("feature", "mtd.list")
    print(nrow(cmp.ft.df.after))
    save(cmp.ft.df.before, cmp.ft.df.after, file = paste0(out.dir, "/", d, "-acrossRank.RData"))
}
stopImplicitCluster()

rm(list = ls())

library(ggplot2)
library(ggupset)

# work.dir <- "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue"
work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/RealData/"

for (d in c("LUADseo", "PRADtcga")) {
    load(paste0(work.dir, "/", d, "-acrossRank.RData"))

    p <- ggplot(cmp.ft.df.before, aes(x = mtd.list)) +
        geom_bar(fill = "#f1a340") +
        scale_x_upset() +
        theme(text = element_text(size = 35, family = "Arial"),
              axis.title.x = element_blank(),
              plot.margin = margin(l = 1.5, t = 0.5, unit = "cm")) +
        theme_combmatrix(combmatrix.label.text = element_text(size = 25, family = "Arial"),
                         combmatrix.panel.line.size = 0.5,
                         combmatrix.label.make_space = F)
    ggsave(paste0(work.dir, "/", d, "-kamratRank-upset-before.svg"), plot = p,
           width = 18, height = 6)

    p <- ggplot(cmp.ft.df.after, aes(x = mtd.list)) +
        geom_bar(fill = "#998ec3") +
        scale_x_upset() +
        theme(text = element_text(size = 35, family = "Arial"),
              axis.title.x = element_blank(),
              plot.margin = margin(l = 1.5, t = 0.5, unit = "cm")) +
        theme_combmatrix(combmatrix.label.text = element_text(size = 25, family = "Arial"),
                         combmatrix.panel.line.size = 0.5,
                         combmatrix.label.make_space = F)
    ggsave(paste0(work.dir, "/", d, "-kamratRank-upset-after.svg"), plot = p,
           width = 18, height = 6)
}
