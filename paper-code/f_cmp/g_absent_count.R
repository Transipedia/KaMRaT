rm(list = ls())

library(stringr)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA"

ci.tab <- NULL
for (w in c("merging-first", "ranking-first")) {
    for (d in c("risk", "relapse")) {
        for (i in 0 : 4) {
            for (m in c("ttest", "snr", "lrc", "nbc")) {
                p <- paste0(work.dir, "/", w, "/", d, "/test", i, "/contig-counts-ridge.", m, ".tsv")
                test.tab <- read.table(p, header = T)
                ci.tab <- rbind(ci.tab, 
                                data.frame("workflow" = w,
                                           "dataset" = d,
                                           "training_set" = i,
                                           "rk.method" = m,
                                           "inserted.nb" = sum(rowSums(test.tab[, -1]) == 0)))                
            }
        }
    }
}

assess.path <- c("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/risk/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/risk/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/relapse/assess.all.tsv")
out.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/cmp_res"

eval.tab <- lapply(assess.path, FUN = function(p) read.table(p, header = T)) %>%
    do.call(what = rbind)

eval.tab$dataset <- str_extract(eval.tab$model_path, pattern = "risk|relapse")
eval.tab$workflow <- str_extract(eval.tab$model_path, pattern = "merging-first|ranking-first")
eval.tab$training_set <- str_extract(eval.tab$model_path, pattern = "train[0-9]+") %>% str_remove(pattern = "train") %>% as.numeric()
eval.tab$rk.method <- str_extract(eval.tab$model_path, pattern = "ttest|snr|lrc|nbc|svm")
eval.tab <- eval.tab[eval.tab$metrics == "balanced_accuracy", ]
eval.tab <- eval.tab[, c("workflow", "dataset", "training_set", "rk.method", "value")]
colnames(eval.tab) <- c("workflow", "dataset", "training_set", "rk.method", "balanced_accuracy")

cmp.df <- merge(x = ci.tab, y = eval.tab, by = c("workflow", "dataset", "training_set", "rk.method"))

ggplot(cmp.df, aes(x = inserted.nb, y = balanced_accuracy, color = workflow)) +
    geom_smooth(method = "lm") +
    geom_point(size = 1) +
    xlab("absent k-mer number") +
    ylab("balanced accuracy") +
    facet_wrap(rk.method ~ ., ncol = 2, nrow = 2, scales = "free") + 
    scale_color_manual(values = c("merging-first" = "royalblue", "ranking-first" = "gold3")) +
    theme(text = element_text(size = 30, family = "Arial"), legend.position = "bottom") +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/cmp_res/bacc_vs_absent.png", width = 12, height = 7)
