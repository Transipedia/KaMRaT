rm(list = ls())

library(stringr)
library(magrittr)
library(ggplot2)

cmp.res <- readLines("/home/haoliang.xue/media/data/PRAD_TCGA_relapse/assess_models.tsv")

cmp.df <- data.frame("run" = str_extract(cmp.res, pattern = "train[0-9]+") %>% str_remove(pattern = "train"),
                     "reduce.mdl" = str_extract(cmp.res, pattern = "ttest|snr|lrc|nbc|svm"),
                     "select.mdl" = str_extract(cmp.res, pattern = "lasso|ridge"),
                     "feature" = str_extract(cmp.res, pattern = "gene|merge-first|rank-first"),
                     "auc" = str_extract(cmp.res, pattern = "0\\.[0-9]+") %>% as.numeric())
ggplot(cmp.df, aes(x = reduce.mdl, y = auc, fill = feature)) +
    geom_hline(yintercept = 0.5, color = "brown") +
    geom_boxplot(width = 0.75, size = 1, outlier.colour = "brown", outlier.shape = "*", outlier.size = 10) +
    facet_grid(. ~ select.mdl) +
    scale_y_continuous(breaks = seq(0.2, 1, 0.1)) +
    scale_fill_manual(values = c("gray", "skyblue", "khaki")) +
    xlab("ranking method") +
    ylab("relapse prediction AUC") +
    theme(text = element_text(size = 26, family = "Arial")) +
    ggsave("/home/haoliang.xue/media/data/PRAD_TCGA_relapse/assess_models.png", width = 12, height = 6)
