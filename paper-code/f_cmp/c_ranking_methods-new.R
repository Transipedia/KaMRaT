rm(list = ls())

library(readxl)
library(stringr)
library(ggplot2)

cmp.df <- read_excel("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/Efficiency.xlsx", sheet = 2)
cmp.df$exec.time <- lapply(cmp.df$CPU_time, 
                           FUN = function(s) {
                               s.l <- str_split(s, pattern = ":") %>% unlist() %>% as.numeric()
                               return((s.l[1] * 3600 + s.l[2] * 60 + s.l[3]) / 3600)
                           }) %>% unlist()
cmp.df$peak_memory <- as.numeric(cmp.df$peak_memory) / 1024 / 1024
cmp.df <- cmp.df[, c("software", "method", "exec.time", "peak_memory")]

cmp.df.long <- pivot_longer(cmp.df, names_to = "group", values_to = "value", cols = c(-software, -method))

p1 <- ggplot(cmp.df.long %>% filter(group == "exec.time"), aes(x = method, y = value + 1, fill = software)) +
    geom_col(color = "black", position = "dodge") +
    scale_y_log10(breaks = c(1, 2, 6, 11, 21, 61, 81), labels = c(0, 1, 5, 10, 20, 60, 80)) +
    scale_fill_manual(values = c("DE-kupl" = "cornsilk3", "KaMRaT" = "lightblue")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())

p2 <- ggplot(cmp.df.long %>% filter(group == "peak_memory"), aes(x = method, y = value + 1, fill = software)) +
    geom_col(color = "black", position = "dodge") +
    scale_y_log10(breaks = c(1, 2, 6, 11, 21, 61, 81), labels = c(0, 1, 5, 10, 20, 60, 80)) +
    scale_fill_manual(values = c("DE-kupl" = "cornsilk3", "KaMRaT" = "lightblue")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())

ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/ranking_methods-new.png", width = 25, height = 16)
