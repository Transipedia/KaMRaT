rm(list = ls())

library(readxl)
library(stringr)
library(ggplot2)

cmp.df <- read_excel("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/Efficiency.xlsx", sheet = 3)
cmp.df$smp.nb <- ifelse(cmp.df$dataset == "risk", yes = "N=204", no = "N=63") %>% factor(levels = c("N=63", "N=204"))
cmp.df$exec.time <- lapply(cmp.df$CPU_time, 
                           FUN = function(s) {
                               s.l <- str_split(s, pattern = ":") %>% unlist() %>% as.numeric()
                               return((s.l[1] * 3600 + s.l[2] * 60 + s.l[3]) / 3600)
                           }) %>% unlist()
cmp.df$memory <- as.numeric(cmp.df$memory) / 1024 / 1024
cmp.df <- cmp.df[, c("smp.nb", "ranking", "exec.time", "memory")]

cmp.df.long <- pivot_longer(cmp.df, names_to = "group", values_to = "value", cols = c(-smp.nb, -ranking))

p1 <- ggplot(cmp.df.long %>% filter(group == "exec.time"), aes(x = smp.nb, y = value + 1, fill = ranking)) +
    geom_col(color = "black", position = "dodge") +
    scale_y_log10(breaks = c(1, 2, 6, 11, 21, 31, 81), labels = c(0, 1, 5, 10, 20, 30, 80)) +
    scale_fill_manual(values = c("oldlace", "plum", "darksalmon", "lightgreen", "cyan")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())

p2 <- ggplot(cmp.df.long %>% filter(group == "memory"), aes(x = smp.nb, y = value + 1, fill = ranking)) +
    geom_col(color = "black", position = "dodge") +
    scale_y_log10(breaks = c(1, 2, 6, 11, 21, 31, 81), labels = c(0, 1, 5, 10, 20, 30, 80)) +
    scale_fill_manual(values = c("oldlace", "plum", "darksalmon", "lightgreen", "cyan")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())

ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/ranking_methods.png", width = 21, height = 19)
