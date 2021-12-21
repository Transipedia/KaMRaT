rm (list = ls())

library(stringr)
library(tidyr)
library(magrittr)
library(foreach)
library(doParallel)

registerDoParallel(makeCluster(6))

get_rc <- function(seq) {
    seq <- as.character(seq)
    seq_rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq_rc)
}

get_cano <- function(seq) {
    seq <- as.character(seq)
    seq_rc <- get_rc(seq)
    if (seq < seq_rc) {
        return(seq)
    }
    else {
        return(seq_rc)
    }
}

mk_kmer_df <- function(ctg, klen, cano = FALSE) {
    ctg <- as.character(ctg)
    ctg.len <- nchar(ctg)
    if (ctg.len < klen) {
        stop("ERROR: contig length < k-mer length")
    }
    kmer.list <- NULL
    for (i in 1 : (ctg.len - klen + 1)) {
        if (!cano) {
            kmer.list <- c(kmer.list, str_sub(ctg, start = i, end = i + klen - 1))        
        }
        else {
            kmer.list <- c(kmer.list, get_cano(str_sub(ctg, start = i, end = i + klen - 1)))
        }
    }
    return(data.frame("kmer.comp" = kmer.list, "contig" = ctg))
}

work.dir <- "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/"
# work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/"
out.dir <- paste0(work.dir, "RealData/")

cmp.df <- NULL
for (d in c("LUADseo", "PRADtcga")) {
    imoka.ft.before <- paste0(work.dir, d, "/imoka_res/CV0/simple/aggregated.kmers.matrix") %>%
        read.table(header = T) %>%
        rownames()
    imoka.ft.before <- imoka.ft.before[-1]
    
    imoka.ft.after <- paste0(work.dir, d, "/imoka_res/CV0/simple/randomforest_models/0_RF.pickle.features") %>%
        read.table(header = F)
    imoka.ft.after <- imoka.ft.after$V1 %>% as.character()
    
    print(all(imoka.ft.before == lapply(imoka.ft.before, get_cano) %>% unlist()))
    print(all(imoka.ft.after == lapply(imoka.ft.after, get_cano) %>% unlist()))
    
    cmp.df <- rbind(cmp.df,
                    foreach (m = c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes"), .combine = rbind) %dopar% {
                        library(magrittr)
                        library(stringr)
                        
                        kamrat.ft.before <- paste0(work.dir, "/", d, "/kamrat_res/CV0/merge-rank/", 
                                                   m, "/top-merged-contig-counts.tsv") %>%
                            read.table(header = T)
                        kamrat.ft.before <- kamrat.ft.before$contig %>% as.character()
                        kamrat.df.before <- lapply(kamrat.ft.before, FUN = function(x) mk_kmer_df(x, 31, T)) %>%
                            do.call(what = rbind)
                        kamrat.df.before.agg <- aggregate(kamrat.df.before$kmer.comp, by = list(kamrat.df.before$contig),
                                                          FUN = function(x) ifelse(length(intersect(as.character(x), imoka.ft.before)) > 0,
                                                                                   yes = intersect(as.character(x), imoka.ft.before)[1],
                                                                                   no = as.character(x)[1]))
                        data.frame("dataset" = d,
                                   "method" = str_remove(m, pattern = "1"),
                                   "group" = "before",
                                   "num.intersect" = length(intersect(imoka.ft.before, 
                                                                      as.character(kamrat.df.before.agg$x))),
                                   "num.imoka" = length(imoka.ft.before),
                                   "num.kamrat" = length(kamrat.ft.before))
                    })
    
    cmp.df <- rbind(cmp.df,
                    foreach (m = c("ttest.padj", "ttest.pi", "snr", "dids", "lr", "bayes"), .combine = rbind) %dopar% {
                        library(magrittr)
                        library(stringr)
                        
                        kamrat.ft.after <- paste0(work.dir, "/", d, "/kamrat_res/CV0/merge-rank/", 
                                                  m, "/randomforest_models/0_RF.pickle.features") %>%
                            read.table(header = F)
                        kamrat.ft.after <- kamrat.ft.after$V1 %>% as.character()
                        kamrat.df.after <- lapply(kamrat.ft.after, FUN = function(x) mk_kmer_df(x, 31, T)) %>%
                            do.call(what = rbind)
                        kamrat.df.after.agg <- aggregate(kamrat.df.after$kmer.comp, by = list(kamrat.df.after$contig),
                                                         FUN = function(x) ifelse(length(intersect(as.character(x), imoka.ft.after)) > 0,
                                                                                  yes = intersect(as.character(x), imoka.ft.after)[1],
                                                                                  no = as.character(x)[1]))
                        data.frame("dataset" = d,
                                   "method" = str_remove(m, pattern = "1"),
                                   "group" = "after",
                                   "num.intersect" = length(intersect(imoka.ft.after, 
                                                                      as.character(kamrat.df.after.agg$x))),
                                   "num.imoka" = length(imoka.ft.after),
                                   "num.kamrat" = length(kamrat.ft.after))
                    })
}
stopImplicitCluster()

cmp.df$ratio.intersect <- cmp.df$num.intersect / (cmp.df$num.imoka + cmp.df$num.kamrat - cmp.df$num.intersect)
write.table(cmp.df, paste0(out.dir, "/cmp_kamrat_imoka.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

rm(list = ls())

library(ggplot2)

out.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/RealData/"
work.path <- paste0(out.dir, "cmp_kamrat_imoka.tsv")

plabeller <- c("LUADseo" = "LUADseo TvsN",
               "PRADtcga" = "PRADtcga RvsNR",
               "before" = "before RF",
               "after" = "after RF")

cmp.df <- read.table(work.path, header = T)
cmp.df$group <- factor(cmp.df$group, levels = c("before", "after"))
p <- ggplot(cmp.df) +
    geom_col(aes(x = method, y = ratio.intersect, fill = group)) +
    geom_text(aes(x = method, y = 0.006, label = paste0(as.character(round(ratio.intersect * 100)), "%")),
              color = "black", size = 7.5) +
    ylab("intersection")+
    facet_grid(group ~ dataset, labeller = as_labeller(plabeller)) +
    scale_y_continuous(breaks = seq(0, 0.2, 0.02), labels = paste0(seq(0, 20, 2), "%"), limits = c(0, 0.15)) +
    scale_fill_manual(values = c("#f1a340", "#998ec3")) +
    theme(text = element_text(size = 25, family = "Arial"),
          legend.position = "none")
ggsave(p, filename = paste0(out.dir, "/cmp_kamrat_imoka.svg"), width = 16, height = 9)
