rm(list = ls())

library(readxl)
library(magrittr)
library(stringr)
library(ggplot2)
library(extrafont)

loadfonts(device = "postscript")

cmp.df <- read_excel("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/Efficiency.xlsx", sheet = 2)
cmp.df$workflow <- factor(cmp.df$workflow, levels = c("ranking-first", "merging-first"))
cmp.df$smp.nb <- ifelse(cmp.df$dataset == "risk", yes = "N=204", no = "N=63") %>% factor(levels = c("N=63", "N=204"))
cmp.df$exec.time <- lapply(cmp.df$CPU_time, 
                           FUN = function(s) {
                               s.l <- str_split(s, pattern = ":") %>% unlist() %>% as.numeric()
                               return((s.l[1] * 3600 + s.l[2] * 60 + s.l[3]) / 3600)
                           }) %>% unlist()
cmp.df$peak_memory <- as.numeric(cmp.df$peak_memory) / 1024 / 1024
cmp.df <- cmp.df[, c("smp.nb", "workflow", "exec.time", "peak_memory")]

cmp.df.long <- pivot_longer(cmp.df, names_to = "group", values_to = "value", cols = c(-smp.nb, -workflow))

p1 <- ggplot(cmp.df.long %>% filter(group == "exec.time"), aes(x = smp.nb, y = value, fill = workflow)) +
    geom_col(color = "black", size = 0.75, position = "dodge") +
    scale_fill_manual(values = c("lightblue", "khaki")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())

p2 <- ggplot(cmp.df.long %>% filter(group == "peak_memory"), aes(x = smp.nb, y = value, fill = workflow)) +
    geom_col(color = "black", size = 0.75, position = "dodge") +
    scale_fill_manual(values = c("lightblue", "khaki")) +
    theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())

ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
    ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/mgfirst_vs_rkfirst.png", width = 12, height = 19)

# to_draw <- group_by(cmp.df.long, smp.nb, workflow, group) %>%
#     mutate(cum_val = cumsum(value))
# 
# p1 <- ggplot(to_draw %>% filter(group == "exec.time"), aes(x = smp.nb, y = cum_val, fill = workflow)) +
#     geom_col(data = . %>% filter(step.num == "2nd"), aes(alpha = step.num), color = "black", size = 0.5, 
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step.num == "1st"), aes(alpha = step.num), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     scale_alpha_manual(values = c("1st" = 0.5, "2nd" = 0), name = "reduction stage") +
#     scale_fill_manual(values = c("navy", "gold3")) +
#     theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank(), axis.text.x = element_blank())
# 
# p2 <- ggplot(to_draw %>% filter(group == "memory"), aes(x = smp.nb, y = cum_val, fill = workflow)) +
#     geom_col(data = . %>% filter(step.num == "2nd"), aes(alpha = step.num), color = "black", size = 0.5, 
#              position = position_dodge(width = 0.9)) +
#     geom_col(data = . %>% filter(step.num == "1st"), aes(alpha = step.num), color = "black", size = 0.5,
#              position = position_dodge(width = 0.9)) +
#     scale_alpha_manual(values = c("1st" = 0.5, "2nd" = 0), name = "reduction stage") +
#     scale_fill_manual(values = c("navy", "gold3")) +
#     theme(text = element_text(size = 65, family = "Arial"), axis.title = element_blank())
# 
# ggarrange(p1, p2, nrow = 2, common.legend = T, align = "v") +
#     ggsave("/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/exec_info/mergingf_vs_rankingf.png", width = 12, height = 19)
