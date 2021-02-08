rm(list = ls())

library(Biostrings)
library(stringr)
library(ggplot2)
library(patchwork)

cmdArgs <- commandArgs(trailingOnly = T)
kamrat.res.dir.prefix <- cmdArgs[1]
spades.res.dir <- cmdArgs[2]
out.dir <- cmdArgs[3]

# kamrat.res.dir.prefix <- "/home/haoliang.xue/media/data/kamrat/paper/f_kamrat_res/a_errfree_randsel/merged-contigs-align-100pct-1-1."
# spades.res.dir <- "/home/haoliang.xue/media/data/kamrat/paper/d_SPAdes_res/a_errfree"
# out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/g_KaMRaT_merge_eval"

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

spades.res <- NULL
for (s in dir(spades.res.dir)) {
    spades.res.x <- read.table(paste0(spades.res.dir, "/", s, "/blastn.alignment.tsv"), header = T)
    spades.res.x$qseqid <- paste0(s, "_", spades.res.x$qseqid)
    spades.res <- rbind(spades.res, spades.res.x)
}
spades.res$align.length <- apply(spades.res[, c("qlen", "length")], MARGIN = 1, FUN = max)
spades.res$rident <- spades.res$nident / spades.res$align.length
spades.res$jacc.idx <- spades.res$nident / (spades.res$qlen + spades.res$slen - spades.res$nident)
spades.res$software <- "SPAdes"

kamrat.res <- NULL
for (m in c("mac_0.32", "pearson_0.21", "spearman_0.23", "none")) {
    kamrat.res.x <- read.table(paste0(kamrat.res.dir.prefix, m, ".tsv"), header = T)
    kamrat.res.x$align.length <- apply(kamrat.res.x[, c("qlen", "length")], MARGIN = 1, FUN = max)
    kamrat.res.x$rident <- kamrat.res.x$nident / kamrat.res.x$align.length
    kamrat.res.x$jacc.idx <- kamrat.res.x$nident / (kamrat.res.x$qlen + kamrat.res.x$slen - kamrat.res.x$nident)
    kamrat.res.x$software <- paste0("KaMRaT-", str_replace(m, pattern = "_", replacement = ":"))
    kamrat.res <- rbind(kamrat.res, kamrat.res.x)
}

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
    scale_fill_manual(values = c("khaki", "azure", "lightblue", "darkseagreen1", "gray")) +
    scale_color_manual(values = c("gold4", "azure4", "royalblue", "darkseagreen4", "dimgray")) +
    ylab("count") +
    theme(text = element_text(size = 22),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

p.box.local <- ggplot(data = to_draw) +
    geom_boxplot(aes(x = len_rg, y = rident, color = software, fill = software), size = 0.75, outlier.size = 1, 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = c("khaki", "azure", "lightblue", "darkseagreen1", "gray")) +
    scale_color_manual(values = c("gold4", "azure4", "royalblue", "darkseagreen4", "dimgray")) +
    ylab("local identical ratio") +
    theme(text = element_text(size = 22),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

p.box.global <- ggplot(data = to_draw) +
    geom_boxplot(aes(x = len_rg, y = jacc.idx, color = software, fill = software), size = 0.75, outlier.size = 1, 
                 position = position_dodge2(width = 0.9, preserve = "single")) +
    scale_fill_manual(values = c("khaki", "azure", "lightblue", "darkseagreen1", "gray")) +
    scale_color_manual(values = c("gold4", "azure4", "royalblue", "darkseagreen4", "dimgray")) +
    xlab("contig length range in group") +
    ylab("jaccard index") +
    theme(text = element_text(size = 22))

p.bar + theme(legend.position = "top") + 
    p.box.local + theme(legend.position = "none") + 
    p.box.global + theme(legend.position = "none") +
    plot_layout(ncol = 1, nrow = 3, heights = c(1, 2, 2)) +
    ggsave(paste0(out.dir, 
                  "/global_local_comparison.", 
                  str_extract(basename(spades.res.dir), pattern = "[a-z]+$"), 
                  str_extract(kamrat.res.dir.prefix, pattern = "-[0-9]+-[0-9]+"), ".png"), 
           width = 21, height = 9)
