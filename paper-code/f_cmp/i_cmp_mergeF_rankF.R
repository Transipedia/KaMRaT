rm (list = ls())

library(Biostrings)
library(magrittr)
library(tidyr)
library(stringr)
library(ggplot2)

getSeqRC <- function(seq)
{
    seq.rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq.rc)
}

retrieveKMers <- function(contig.list, k.len, canonical) {
    tmp.df <- lapply(contig.list,
                     FUN = function(ctg) {
                         seq.len <- nchar(ctg)
                         kmer.vect <- NULL
                         for (start_pos in 1 : (seq.len - k.len + 1)) {
                             kmer.seq <- substr(unname(ctg), start_pos, start_pos + k.len - 1)
                             if (canonical) {
                                 kmer.seq.rc <- getSeqRC(kmer.seq)
                                 if (kmer.seq.rc < kmer.seq) {
                                     kmer.seq <- kmer.seq.rc
                                 }
                             }
                             kmer.vect <- c(kmer.vect, kmer.seq)
                         }
                         return(data.frame("contig" = ctg, "kmer" = kmer.vect))
                     }) %>%
        do.call(what = rbind)
    return(tmp.df)
}

rankf.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/"
mergef.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/"
cmp.summary <- NULL
for (d in c("relapse", "risk")) {
    for (i in 0 : 4) {
        for (m in c("ttest", "snr", "lrc", "nbc")) {
            print(paste0(d, "-", i, "-", m))

            rankf.path <- paste0(rankf.prefix, d, "/train", i, "/", m, "/top5000-contig-counts.tsv")
            rankf.res <- read.table(rankf.path, header = T)$contig %>% as.character()
            rankf.ctg.kmer.df <- retrieveKMers(rankf.res, 31, T)
            
            mergef.path <- paste0(mergef.prefix, d, "/train", i, "/", m, "/top5000-contig-counts.tsv")
            mergef.res <- read.table(mergef.path, header = T)$contig %>% as.character()
            mergef.ctg.kmer.df <- retrieveKMers(mergef.res, 31, T)
            
            shared.ctg <- merge(x = rankf.ctg.kmer.df, y = mergef.ctg.kmer.df, by = "kmer")
            shared.nb <- min(length(unique(shared.ctg$contig.x)), length(unique(shared.ctg$contig.y)))
            cmp.summary <- rbind(cmp.summary,
                                 data.frame("dataset" = d,
                                            "train_set" = paste0("train", i),
                                            "method" = m,
                                            "shared.nb" = shared.nb))
        }
        
    }
}

cmp.summary.mean <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = mean)
colnames(cmp.summary.mean) <- c("dataset", "method", "shared.mean")
cmp.summary.max <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = max)
colnames(cmp.summary.max) <- c("dataset", "method", "shared.max")
cmp.summary.min <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = min)
colnames(cmp.summary.min) <- c("dataset", "method", "shared.min")

cmp.summary.agg <- merge(cmp.summary.mean, cmp.summary.max) %>%
    merge(cmp.summary.min)


ggplot(cmp.summary.agg) +
    geom_col(aes(x = dataset, y = shared.mean, fill = "method"), color = "black", size = 1) +
    geom_errorbar(aes(x = dataset, ymin = shared.min, ymax = shared.max), color = "black", size = 1, width = 0.5) +
    theme(text = element_text(size = 35, family = "Arial")) +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/new_cmp_res/mergef_vs_rankf.png", 
           width = 8, height = 11)
