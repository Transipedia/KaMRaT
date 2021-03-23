rm(list = ls())

library(stringr)
library(rjson)
library(ggplot2)

assess.path <- c("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/risk/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse-new/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/risk/assess.all.tsv",
                 "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/relapse-new/assess.all.tsv")
out.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/cmp_res"

eval.tab <- lapply(assess.path, FUN = function(p) read.table(p, header = T)) %>%
    do.call(what = rbind)

eval.tab$dataset <- ifelse(str_detect(eval.tab$model_path, pattern = "risk"), yes = "risk (N=204)", no = "relapse (N=63)")
eval.tab$workflow <- str_extract(eval.tab$model_path, pattern = "merging-first|ranking-first")
eval.tab$train.set <- str_extract(eval.tab$model_path, pattern = "train[0-9]+")
eval.tab$rk.method <- str_extract(eval.tab$model_path, pattern = "ttest|snr|lrc|nbc|svm")
eval.tab$rk.method <- factor(eval.tab$rk.method, levels = c("ttest", "snr", "lrc", "nbc", "svm"))

eval.tab <- eval.tab[eval.tab$metrics == "balanced_accuracy", ]
eval.tab <- eval.tab[eval.tab$rk.method != "svm", ]
eval.tab$software <- "KaMRaT"
eval.tab <- eval.tab[, c("train.set", "rk.method", "value", "workflow", "dataset", "software")]

imoka.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA"
imoka.bacc <- lapply(0 : 4,
                     FUN = function(i) {
                         bacc.risk <- fromJSON(file = paste0(imoka.prefix, "/risk/test", i, "/prediction_test.txt"))$balanced_accuracy
                         bacc.relapse <- fromJSON(file = paste0(imoka.prefix, "/relapse/test", i, "/prediction_test.txt"))$balanced_accuracy
                         return(data.frame("train.set" = i, 
                                           "value" = c(bacc.risk, bacc.relapse), 
                                           "dataset" = c("risk (N=204)", "relapse (N=63)")))
                     }) %>%
    do.call(what = rbind)
imoka.bacc$software <- "iMOKA"
imoka.bacc$rk.method <- "iMOKA"
imoka.bacc$workflow <- "iMOKA"

cmp.df <- rbind(eval.tab, imoka.bacc)
cmp.df$software <- factor(cmp.df$software, levels = c("iMOKA", "KaMRaT"))
cmp.df$workflow <- factor(cmp.df$workflow, levels = c("iMOKA", "merging-first", "ranking-first"))

ggplot(cmp.df, aes(x = rk.method, y = value, fill = workflow)) +
    geom_hline(yintercept = 0.5, color = "dimgray", linetype = "dashed") +
    geom_boxplot() +
    # geom_dotplot(binaxis = 'y', binwidth = 0.01, dotsize = 3, stackdir = "center", position = position_dodge(0.35)) +
    facet_grid(dataset ~ software, scales = "free_x", space = "free_x") +
    scale_y_continuous(breaks = seq(0.4, 1, 0.1)) +
    scale_fill_manual(values = c("merging-first" = "lightblue", "ranking-first" = "khaki", "iMOKA" = "gray")) +
    xlab("ranking methods") +
    ylab("balanced accuracy") +
    theme(text = element_text(size = 30, family = "Arial"), legend.position = "top", axis.title.x = element_blank()) +
    ggsave(paste0(out.dir, "/balanced_accuracy-1-5.png"), width = 16, height = 7)
