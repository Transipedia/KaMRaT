rm (list = ls())

library(readxl)
library(stringr)
library(ggplot2)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/RealData/"

perf.tab <- read_excel(paste0(work.dir, "KaMRaT vs iMOKA newest.xlsx"), sheet = 2)
perf.tab$median <- apply(perf.tab[, str_detect(colnames(perf.tab), pattern = "train")], MARGIN = 1, 
                         FUN = function(x) median(as.numeric(x)))
perf.tab$min <- apply(perf.tab[, str_detect(colnames(perf.tab), pattern = "train")], MARGIN = 1, 
                      FUN = function(x) min(as.numeric(x)))
perf.tab$max <- apply(perf.tab[, str_detect(colnames(perf.tab), pattern = "train")], MARGIN = 1, 
                      FUN = function(x) max(as.numeric(x)))
perf.tab$method <- factor(perf.tab$method, 
                          levels = c("iMOKA", "bayes", "dids", "lr", "snr", "ttest.padj", "ttest.pi"))

p <- ggplot(perf.tab) +
    geom_col(aes(x = method, y = median, fill = workflow), position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(x = method, ymin = min, ymax = max, group = workflow), width = 0.5, color = "black", 
                  position = position_dodge(width = 0.9)) +
    facet_grid(dataset ~ .) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_fill_manual(values = c("iMOKA-default" = "#5e3c99", "iMOKA-simple" = "#998ec3", 
                                 "KaMRaT-merge-rank" = "#fdb863", "KaMRaT-rank-merge" = "#e66101")) +
    ylab("balanced accuracy") +
    theme(text = element_text(size = 25, family = "Arial"), legend.position = "none",
          axis.title.x = element_blank())
ggsave(paste0(work.dir, "/kamrat-imoka-perf.svg"), plot = p, width = 9, height = 8)

effi.tab <- read_excel(paste0(work.dir, "KaMRaT vs iMOKA newest.xlsx"), sheet = 1)
effi.tab$exec.time <- lapply(effi.tab$CPU_time, 
                             FUN = function(s) {
                                 s.l <- str_split(s, pattern = ":") %>% unlist() %>% as.numeric()
                                 return((s.l[1] * 3600 + s.l[2] * 60 + s.l[3]) / 3600)
                             }) %>% unlist()
effi.tab$memory <- as.numeric(effi.tab$memory) / 1024 / 1024
effi.tab$method <- factor(effi.tab$method, 
                          levels = c("iMOKA", "bayes", "dids", "lr", "snr", "ttest.padj", "ttest.pi"))

effi.tab.time <- aggregate(effi.tab$exec.time,
                           by = list(effi.tab$dataset, effi.tab$workflow, effi.tab$method, effi.tab$set),
                           FUN = sum)
colnames(effi.tab.time) <- c("dataset", "workflow", "method", "set", "exec.time")
effi.tab.time2 <- aggregate(effi.tab.time$exec.time, 
                            by = list(effi.tab.time$dataset, effi.tab.time$workflow, effi.tab.time$method), 
                            FUN = mean)
colnames(effi.tab.time2) <- c("dataset", "workflow", "method", "exec.time")

p <- ggplot(effi.tab.time2) +
    geom_col(aes(x = method, y = exec.time, fill = workflow), color = "black", position = position_dodge()) +
    facet_grid(dataset ~ ., scales = "free") +
    scale_fill_manual(values = c("iMOKA-simple" = "#998ec3", 
                                 "KaMRaT-merge-rank" = "#fdb863", "KaMRaT-rank-merge" = "#e66101")) +
    ylab("CPU time (hours)") +
    theme(text = element_text(size = 25, family = "Arial"), legend.position = "none",
          axis.title.x = element_blank())
ggsave(paste0(work.dir, "/kamrat-imoka-cput.svg"), plot = p, width = 9, height = 8)

effi.tab.mem <- aggregate(effi.tab$memory, 
                          by = list(effi.tab$dataset, effi.tab$workflow, effi.tab$method), 
                          FUN = max)
colnames(effi.tab.mem) <- c("dataset", "workflow", "method", "peak.mem")
p <- ggplot(effi.tab.mem) +
    geom_col(aes(x = method, y = peak.mem, fill = workflow), color = "black", position = position_dodge()) +
    facet_grid(dataset ~ .) +
    scale_fill_manual(values = c("iMOKA-simple" = "#998ec3", 
                                 "KaMRaT-merge-rank" = "#fdb863", "KaMRaT-rank-merge" = "#e66101")) +
    ylab("peak RAM (GB)") +
    theme(text = element_text(size = 25, family = "Arial"), legend.position = "none",
          axis.title.x = element_blank())
ggsave(paste0(work.dir, "/kamrat-imoka-mem.svg"), plot = p, width = 9, height = 8)
