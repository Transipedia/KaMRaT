rm(list = ls())

library(readxl)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(extrafont)

loadfonts(device = "postscript")

cmp.df <- read_excel("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/Efficiency.xlsx", sheet = 1)
cmp.df$exec.time <- lapply(cmp.df$CPU_time, 
                                  FUN = function(s) {
                                      s.l <- str_split(s, pattern = ":") %>% unlist() %>% as.numeric()
                                      return((s.l[1] * 3600 + s.l[2] * 60 + s.l[3]) / 3600)
                                  }) %>% unlist()
cmp.df$peak_memory <- as.numeric(cmp.df$peak_memory) / 1024 / 1024
cmp.df <- cmp.df[, c("software", "workflow", "exec.time", "peak_memory")]

cmp.df.long <- pivot_longer(cmp.df, names_to = "group", values_to = "value", cols = c(-software, -workflow))

p1 <- ggplot(cmp.df.long %>% filter(group == "exec.time"), aes(x = software, y = value, fill = workflow)) +
    geom_col(color = "black", size = 0.75, position = "dodge") +
    scale_fill_manual(values = c("azure3", "lightblue", "khaki")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())

p2 <- ggplot(cmp.df.long %>% filter(group == "peak_memory"), aes(x = software, y = value, fill = workflow)) +
    geom_col(color = "black", size = 0.75, position = "dodge") +
    scale_fill_manual(values = c("azure3", "lightblue", "khaki")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())

ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/imoka_vs_kamrat-new.png", width = 12, height = 16)

# to_draw <- group_by(cmp.df.long, smp.nb, software, group) %>%
#     mutate(cum_val = cumsum(value))
# 
# p1 <- ggplot(to_draw %>% filter(group == "exec.time"), aes(x = smp.nb, y = cum_val, fill = software)) +
#     geom_col(data = . %>% filter(step == "2-extend"), aes(alpha = step), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step == "1-reduce"), aes(alpha = step), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step == "0-joinMat"), aes(alpha = step), color = "black", size = 0.5, 
#              position = position_dodge(width = 0.9)) +
#     scale_alpha_manual(values = c("2-extend" = 0, "1-reduce" = 0.5, "0-joinMat" = 1)) +
#     scale_fill_manual(values = c("darkred", "navy")) +
#     theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())
# 
# p2 <- ggplot(to_draw %>% filter(group == "memory"), aes(x = smp.nb, y = cum_val, fill = software)) +
#     geom_col(data = . %>% filter(step == "2-extend"), aes(alpha = step), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step == "1-reduce"), aes(alpha = step), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step == "0-joinMat"), aes(alpha = step), color = "black", size = 0.5, 
#              position = position_dodge(width = 0.9)) +
#     scale_alpha_manual(values = c("2-extend" = 0, "1-reduce" = 0.5, "0-joinMat" = 1)) +
#     scale_fill_manual(values = c("darkred", "navy")) +
#     theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())
# 
# ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
#     ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/imoka_vs_kamrat.png", width = 12, height = 19)
