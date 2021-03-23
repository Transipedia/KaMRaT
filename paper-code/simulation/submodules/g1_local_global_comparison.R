rm(list = ls())

library(Biostrings)
library(stringr)
library(ggplot2)
library(patchwork)

cmdArgs <- commandArgs(trailingOnly = T)
kamrat.res.dir.prefix <- cmdArgs[1]
spades.res.dir <- cmdArgs[2]
out.dir <- cmdArgs[3]

kamrat.res.dir.prefix <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/simulation/f_kamrat_res/a_randsel/errfree/merged-contigs-align-100pct-1-1."
spades.res.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/simulation/d_SPAdes_res/a_errfree"
out.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/simulation/g_KaMRaT_merge_eval"

get_quantile <- function(val, val.quant.vect) {
    bin.left <- val.quant.vect[-length(val.quant.vect)]
    bin.right <- val.quant.vect[-1]
    if (val == bin.left[1]) {
        return(1)
    }
    else {
        return(min(which(val > bin.left & val <= bin.right)))
    }
}

spades.res <- read.table(paste0(spades.res.dir, "/all_samples/blastn.alignment.tsv"), header = T)
spades.res$qseqid <- paste0("all_samples", "_", spades.res$qseqid)
spades.res$align.length <- apply(spades.res[, c("qlen", "length")], MARGIN = 1, FUN = max)
spades.res$rident <- spades.res$nident / spades.res$align.length
spades.res$jacc.idx <- spades.res$nident / (spades.res$qlen + spades.res$slen - spades.res$nident)
spades.res$software <- "SPAdes"

kamrat.res <- read.table(paste0(kamrat.res.dir.prefix, "none.tsv"), header = T)
kamrat.res$align.length <- apply(kamrat.res[, c("qlen", "length")], MARGIN = 1, FUN = max)
kamrat.res$rident <- kamrat.res$nident / kamrat.res$align.length
kamrat.res$jacc.idx <- kamrat.res$nident / (kamrat.res$qlen + kamrat.res$slen - kamrat.res$nident)
kamrat.res$software <- "KaMRaT-none"

to_draw <- rbind(spades.res[, c("qlen", "rident", "jacc.idx", "software")], 
                 kamrat.res[, c("qlen", "rident", "jacc.idx", "software")])
to_draw.len.qt <- unname(quantile(to_draw$qlen, probs = seq(from = 0, to = 1, by = 0.1)))
to_draw$n_qt <- unlist(lapply(to_draw$qlen, FUN = function(x) get_quantile(x, to_draw.len.qt)))

to_draw$len_rg <- NA
for (i_qt in 1 : 10) {
    len.list <- to_draw$qlen[to_draw$n_qt == i_qt]
    to_draw$len_rg[to_draw$n_qt == i_qt] <- paste0(min(len.list), "-", max(len.list))
}
to_draw$len_rg <- reorder(x = to_draw$len_rg, X = to_draw$n_qt)

to_draw.sum <- aggregate(x = to_draw$len_rg, by = list(to_draw$len_rg, to_draw$software), FUN = function(x) length(x))
names(to_draw.sum) <- c("len_rg", "software", "count")

p.bar <- ggplot(data = to_draw.sum) +
    geom_col(aes(x = len_rg, y = count, fill = software, color = software, group = software), 
             size = 0.75, position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = c("azure", "gray")) +
    scale_color_manual(values = c("azure4", "dimgray")) +
    ylab("count") +
    theme(text = element_text(size = 26, family = "Arial"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

p.box.local <- ggplot(data = to_draw) +
    geom_boxplot(aes(x = len_rg, y = rident, color = software, fill = software), size = 0.75, outlier.size = 1, 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = c("azure", "gray")) +
    scale_color_manual(values = c("azure4", "dimgray")) +
    ylab("local identical ratio") +
    theme(text = element_text(size = 26, family = "Arial"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

p.box.global <- ggplot(data = to_draw) +
    geom_boxplot(aes(x = len_rg, y = jacc.idx, color = software, fill = software), size = 0.75, outlier.size = 1, 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = c("azure", "gray")) +
    scale_color_manual(values = c("azure4", "dimgray")) +
    xlab("contig length range in group") +
    ylab("jaccard index") +
    theme(text = element_text(size = 26, family = "Arial"))

p.bar + theme(legend.position = "top") + 
    p.box.local + theme(legend.position = "none") + 
    p.box.global + theme(legend.position = "none") +
    plot_layout(ncol = 1, nrow = 3, heights = c(1, 2, 2)) +
    ggsave(paste0(out.dir, 
                  "/global_local_comparison.", 
                  str_extract(basename(spades.res.dir), pattern = "[a-z]+$"), 
                  str_extract(kamrat.res.dir.prefix, pattern = "-[0-9]+-[0-9]+"), ".png"), 
           width = 16, height = 10)
