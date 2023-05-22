rm (list = ls())

library(stringr)
library(ggplot2)
library(ggfortify)
library(patchwork)

# Sample number vs k-mer number
cmp.tab.kmer <- read.csv("../../../../revision/merge-on-real/CPUtime_and_peakRAM/kmerNumber_and_sampleNumber.csv",
                         header = TRUE)
cmp.tab.kmer$IndexSize.gb <- (str_remove(cmp.tab.kmer$Index.Size, "kb") %>% as.integer()) / 1000 / 1000

plt.kmernum <- ggplot(data = cmp.tab.kmer) +
    geom_point(aes(x = Sample.Number, y = k.mer.Number), size = 3) +
    geom_line(aes(x = Sample.Number, y = k.mer.Number), size = 1) +
    scale_x_continuous(breaks = c(10, 20, 50, 80, 140)) +
    scale_y_continuous(breaks = seq(0, 5e8, 5e7)) +
    xlab("sample number") +
    ylab("filtered k-mer number") +
    theme_bw() +
    theme(text = element_text(size = 15), panel.grid = element_blank())
ggsave(plt.kmernum,
       filename = "../../../../revision/merge-on-real/CPUtime_and_peakRAM/kmerNumber_and_sampleNumber.pdf",
       width = 6, height = 9)

plt.indexsize <- ggplot(data = cmp.tab.kmer) +
    geom_point(aes(x = Sample.Number, y = IndexSize.gb), colour = "#5e3c99", size = 3) +
    geom_line(aes(x = Sample.Number, y = IndexSize.gb), colour = "#5e3c99", size = 1) +
    scale_x_continuous(breaks = c(10, 20, 50, 80, 140)) +
    scale_y_continuous(breaks = seq(0, 350, 50)) +
    xlab("sample number") +
    ylab("kamrat index size") +
    theme_bw() +
    theme(text = element_text(size = 15), panel.grid = element_blank())
ggsave(plt.indexsize,
       filename = "../../../../revision/merge-on-real/CPUtime_and_peakRAM/indexSize_and_sampleNumber.pdf",
       width = 6, height = 9)


# KaMRaT performance
cmp.tab.kamrat <- read.csv("../../../../revision/merge-on-real/CPUtime_and_peakRAM/KaMRaT_CPUtime_and_peakRAM.csv",
                           header = TRUE)
cmp.tab.kamrat$peakRAM.gb <- (str_remove(cmp.tab.kamrat$Peak.RAM, "kb") %>% as.integer()) / 1000 / 1000
cmp.tab.kamrat$CPUtime.h <- lapply(cmp.tab.kamrat$CPU.Time,
                                   FUN = function(x) {
                                       hms <- (str_split(x, ":")[[1]] %>% as.integer())
                                       return((hms[1] * 3600 + hms[2] * 60 + hms[3])/3600)
                                   }) %>% unlist()
## Scaling curves
cmp.tab.scal <- cmp.tab.kamrat[cmp.tab.kamrat$Mode %in% c("with NF", "pearson", "snr"), ]
cmp.tab.scal$CmdMode <- paste(cmp.tab.scal$Command, cmp.tab.scal$Mode, sep = ": ")
plt.cput <- ggplot(data = cmp.tab.scal) +
    # geom_smooth(aes(x = Sample.Number, y = CPUtime.h, colour = CmdMode, fill = CmdMode),
    #             formula = y ~ x, size = 1, method = "glm") +
    geom_line(aes(x = Sample.Number, y = CPUtime.h, colour = CmdMode), size = 1) +
    geom_point(aes(x = Sample.Number, y = CPUtime.h, colour = CmdMode), size = 3) +
    scale_colour_manual(values = c("index: with NF" = "#5e3c99",
                                   "rank: snr" = "#e66101",
                                   "merge: pearson" = "#fdb863")) +
    scale_fill_manual(values = c("index: with NF" = "#5e3c99",
                                 "rank: snr" = "#e66101",
                                 "merge: pearson" = "#fdb863")) +
    scale_x_continuous(breaks = c(10, 20, 50, 80, 140, 170, 232)) +
    scale_y_continuous(breaks = seq(0, 50, 5)) +
    xlab(NULL) +
    ylab("CPU time (h)") +
    theme_bw() +
    theme(text = element_text(size = 15), axis.text.x = element_blank(),
          panel.grid = element_blank())
plt.ram <- ggplot(data = cmp.tab.scal) +
    # geom_smooth(aes(x = Sample.Number, y = peakRAM.gb, colour = CmdMode, fill = CmdMode),
    #             formula = y ~ x, size = 1, method = "glm") +
    geom_line(aes(x = Sample.Number, y = peakRAM.gb, colour = CmdMode), size = 1) +
    geom_point(aes(x = Sample.Number, y = peakRAM.gb, colour = CmdMode), size = 3) +
    scale_colour_manual(values = c("index: with NF" = "#5e3c99",
                                   "rank: snr" = "#e66101",
                                   "merge: pearson" = "#fdb863")) +
    scale_fill_manual(values = c("index: with NF" = "#5e3c99",
                                 "rank: snr" = "#e66101",
                                 "merge: pearson" = "#fdb863")) +
    scale_x_continuous(breaks = c(10, 20, 50, 80, 140, 170, 232)) +
    scale_y_continuous(breaks = seq(0, 300, 50)) +
    xlab("sample number") +
    ylab("peak RAM (gb)") +
    theme_bw() +
    theme(text = element_text(size = 15), panel.grid = element_blank())
plt <- ((plt.cput / plt.ram) & theme(legend.position = "right")) + 
    plot_layout(guides = "collect")
ggsave(plt, width = 9, height = 9,
       filename = "../../../../revision/merge-on-real/CPUtime_and_peakRAM/KaMRaT_CPUtime_and_peakRAM_scal.pdf")

## Bar plots in situation 140 samples
cmp.tab.140smp <- cmp.tab.kamrat[cmp.tab.kamrat$Sample.Number == 140, ]
plt.lst <- lapply(c("index", "merge", "rank"),
                  FUN = function(cmd) {
                      cmp.tab2plot <- cmp.tab.140smp[cmp.tab.140smp$Command == cmd, ]
                      print(cmd)
                      plt.cputime <- ggplot(data = cmp.tab2plot) +
                          geom_col(aes(x = Mode, y = CPUtime.h, fill = Command)) +
                          scale_fill_manual(values = c("index" = "#5e3c99",
                                                       "rank" = "#e66101",
                                                       "merge" = "#fdb863")) +
                          scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
                          xlab(NULL) +
                          ylab("CPU time (h)") +
                          theme_bw() +
                          theme(text = element_text(size = 15), axis.text.x = element_blank(),
                                panel.grid = element_blank(), legend.position = "none")
                      plt.peakRAM <- ggplot(data = cmp.tab2plot) +
                          geom_col(aes(x = Mode, y = peakRAM.gb, fill = Command)) +
                          scale_fill_manual(values = c("index" = "#5e3c99",
                                                       "rank" = "#e66101",
                                                       "merge" = "#fdb863")) +
                          scale_y_continuous(breaks = seq(0, 300, 50), limits = c(0, 300)) +
                          xlab(NULL) +
                          ylab("peak RAM (gb)") +
                          theme_bw() +
                          theme(text = element_text(size = 15), panel.grid = element_blank(),
                                legend.position = "none")
                      if (cmd == "index") {
                          return(plt.cputime / plt.peakRAM)
                      } else {
                          return(plt.cputime / plt.peakRAM & 
                                     theme(axis.title.y = element_blank(),
                                           axis.text.y = element_blank()))
                      }
                  })
plt <- (plt.lst[[1]] | plt.lst[[2]] | plt.lst[[3]]) + plot_layout(width = c(2.5, 4, 6))
ggsave(plt, width = 12, height = 9,
       filename = "../../../../revision/merge-on-real/CPUtime_and_peakRAM/KaMRaT_CPUtime_and_peakRAM_140smp.pdf")
