rm(list = ls())

library(stringr)
library(tidyr)
library(magrittr)
library(ggplot2)

nb.summary <- read.table("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/merging-ratio.tsv", header = F)
colnames(nb.summary) <- c("number", "file")
nb.summary$number <- nb.summary$number - 1
nb.summary$type <- ifelse(str_detect(nb.summary$file, pattern = "CVMatrices"), yes = "kmer", no = "contig")
nb.summary$dataset <- str_extract(nb.summary$file, pattern = "risk|relapse")
nb.summary$train_set <- str_extract(nb.summary$file, pattern = "train[0-9]+")
nb.summary <- nb.summary[, c("dataset", "train_set", "type", "number")]

nb.summary.wide <- pivot_wider(nb.summary, id_cols = c("dataset", "train_set"), names_from = "type", values_from = "number")
nb.summary.wide$reduc.ratio <- nb.summary.wide$kmer / nb.summary.wide$contig

nb.summary.mean <- aggregate(nb.summary.wide$reduc.ratio, by= list(nb.summary.wide$dataset), FUN = mean)
colnames(nb.summary.mean) <- c("dataset", "mean")
nb.summary.min <- aggregate(nb.summary.wide$reduc.ratio, by= list(nb.summary.wide$dataset), FUN = min)
colnames(nb.summary.min) <- c("dataset", "min")
nb.summary.max <- aggregate(nb.summary.wide$reduc.ratio, by= list(nb.summary.wide$dataset), FUN = max)
colnames(nb.summary.max) <- c("dataset", "max")
nb.summary.agg <- merge(nb.summary.mean, nb.summary.max) %>%
    merge(nb.summary.min)

ggplot(nb.summary.agg) +
    geom_col(aes(x = dataset, y = mean), color = "black", size = 1, fill = "gray") +
    geom_errorbar(aes(x = dataset, ymin = min, ymax = max), size = 1, width = 0.5) +
    ylab("reduction ratio") +
    theme(text = element_text(size = 35, family = "Arial")) +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/cmp_res/reduction_ratio.png", width = 8, height = 11)
