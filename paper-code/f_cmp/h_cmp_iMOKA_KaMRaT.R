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
    tmp.df <- NULL
    for (ctg in contig.list) {
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
        tmp.df <- rbind(tmp.df, 
                        data.frame("contig" = ctg, "kmer" = kmer.vect))
    }
    return(tmp.df)
}

imoka.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA/"
kamrat.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/"
cmp.summary <- NULL
for (d in c("relapse", "risk")) {
    for (i in 0 : 4) {
        print(paste0(d, "-", i))
        imoka.path <- paste0(imoka.prefix, d, "/train", i, "/randomforest_models/0_RF.pickle.features")
        imoka.res <- readLines(imoka.path)
        imoka.res.rc <- lapply(imoka.res, function(s) getSeqRC(s)) %>% unlist()
        kamrat.path <- paste0(kamrat.prefix, d, "/train", i, "/nbc/top5000-contig-counts.tsv")
        kamrat.res <- read.table(kamrat.path, header = T)$contig %>% as.character()
        kamrat.ctg.kmer.df <- retrieveKMers(kamrat.res, 31, T)
        shared.nb <- length(unique(kamrat.ctg.kmer.df[kamrat.ctg.kmer.df$kmer %in% c(imoka.res, imoka.res.rc), "contig"]))
        cmp.summary <- rbind(cmp.summary,
                             data.frame("dataset" = d,
                                        "train_set" = paste0("train", i),
                                        "shared.nb" = shared.nb))
    }
}

cmp.summary.mean <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = mean)
colnames(cmp.summary.mean) <- c("dataset", "shared.mean")
cmp.summary.max <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = max)
colnames(cmp.summary.max) <- c("dataset", "shared.max")
cmp.summary.min <- aggregate(cmp.summary$shared.nb, by = list(cmp.summary$dataset), FUN = min)
colnames(cmp.summary.min) <- c("dataset", "shared.min")

cmp.summary.agg <- merge(cmp.summary.mean, cmp.summary.max) %>%
    merge(cmp.summary.min)


ggplot(cmp.summary.agg) +
    geom_col(aes(x = dataset, y = shared.mean), color = "black", fill = "gray", size = 1) +
    geom_errorbar(aes(x = dataset, ymin = shared.min, ymax = shared.max), color = "black", size = 1, width = 0.5) +
    theme(text = element_text(size = 35, family = "Arial")) +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/new_cmp_res/kamrat_vs_imoka.png", 
           width = 8, height = 11)
