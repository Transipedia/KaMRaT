rm(list = ls())

library(Biostrings)
library(stringr)
library(ggplot2)

cmdArg <- commandArgs(trailingOnly = T)
spades.dir <- cmdArg[1]
raw.dir.prefix <- cmdArg[2]
kamrat.dir.prefix <- cmdArg[3]
out.dir <- cmdArg[4]

spades.dir <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/d_SPAdes_res/a_errfree"
raw.dir.prefix <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/e_jellyfish_joinCount_res/a_errfree/31/raw-counts-"
kamrat.dir.prefix <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/f_kamrat_res"
out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/simulation/g_KaMRaT_merge_eval"

set.seed(91400)

calc_wrong_ratio <- function(align.tsv.path, all.fa.path) {
    aligntab <- read.table(align.tsv.path, header = T)
    aligntab$is.correct <- (aligntab$qlen == aligntab$qend & 
                                         aligntab$qstart == 1 &
                                         aligntab$qlen == aligntab$nident &
                                         aligntab$pident == 100)
    fa.all <- readDNAStringSet(all.fa.path)
    fa.all <- fa.all[width(fa.all) > 31]
    fa.wrong <- fa.all[!(names(fa.all) %in% aligntab$qseqid[aligntab$is.correct])]
    return(length(fa.wrong) / length(fa.all))
}

## First Comparison: SPAdes vs KaMRaT-merge ##
spades.res.prefix <- paste0(spades.dir, "/all_samples")
spades.res <- calc_wrong_ratio(paste0(spades.res.prefix, "/blastn.alignment.tsv"),
                               paste0(spades.res.prefix, "/transcripts.fasta"))
kamrat.none.prefix <- paste0(kamrat.dir.prefix, "/a_errfree_randsel")
kamrat.res <- calc_wrong_ratio(paste0(kamrat.none.prefix, "/merged-contigs-align-100pct-1-1.none.tsv"),
                               paste0(kamrat.none.prefix, "/merged-contigs-100pct-1-1.none.fa"))
out.file <- file(paste0(out.dir, "/wrong_contig_ratio_KaMRaT_vs_SPAdes.txt"))
writeLines(c(paste0("SPAdes wrong assembly ratio: ", round(spades.res, 5) * 100, "%"),
             paste0("KaMRaT-none wrong extension ratio: ", round(kamrat.res, 5) * 100, "%")), out.file)
close(out.file)

## Figure 1: Comparison with variant k ##
kamrat.dir <- paste0(kamrat.dir.prefix, "/b_errfree_interval")
extens.qual <- NULL
for (k in c(31, 25, 21, 19, 15)) {
    for (m in seq(k - 1, floor(k / 2), by = -1)) {
        for (x in c("none", "mac_0.32", "pearson_0.21", "spearman_0.23")) {
            print(paste0(k, "-", m, ": ", x))
            aligntab.path <- paste0(kamrat.dir, "/", k, "-", m, "/merged-contigs-align-100pct-1-1.", x, ".tsv")
            fa.path <- paste0(kamrat.dir, "/", k, "-", m, "/merged-contigs-100pct-1-1.", x, ".fa")
            extens.qual <- rbind(extens.qual,
                                 data.frame("k.len" = k,
                                            "min.overlap" = m,
                                            "method" = str_replace(x, pattern = "_", replacement = ":"),
                                            "wrong.ratio" = calc_wrong_ratio(aligntab.path, fa.path)))
        }
    }
}

ggplot(extens.qual) +
    geom_line(aes(x = min.overlap, y = wrong.ratio, color = method), size = 1.5) +
    facet_grid(paste0("k=", k.len) ~ ., scales = "fixed") +
    scale_color_manual(values = c("azure4", "gold3", "royalblue", "darkseagreen3")) +
    scale_x_reverse() +
    xlab("min overlap length") +
    ylab("wrong extension ratio") +
    theme(text = element_text(size = 26),
	  legend.position = "top") +
    ggsave(paste0(out.dir, "/wrong_extension_ratio_with_varklen.png"), width = 16, height = 10)

rm(extens.qual, k, kamrat.dir, m, aligntab.path, fa.path, x)

## Figure 2: Comparison with different kept ratios ##
nb.kmer.all <- nrow(read.table(paste0(raw.dir.prefix, "100pct-1-1.tsv"), header = T))
extens.qual <- NULL
for (d in c("errfree", "witherr")) {
    kamrat.dir <- paste0(kamrat.dir.prefix, "/a_randsel/", d)
    for (r in seq(100, 10, by = -10)) {
        for (x in c("none", "mac_0.32", "pearson_0.21", "spearman_0.23")) {
            print(paste0(d, "-", r, "-", x))
            raw.kmer.nb <- nrow(read.table(paste0(raw.dir.prefix, r, "pct-1-1.tsv"), header = T))
            if (d == "errfree") {
                a <- 1
            } else {
                a <- 3
            }
            aligntab.path <- paste0(kamrat.dir, "/merged-contigs-align-", r, "pct-1-", a, ".", x, ".tsv")
            fa.path <- paste0(kamrat.dir, "/merged-contigs-", r, "pct-1-", a, ".", x, ".fa")
            extens.qual <- rbind(extens.qual,
                                 data.frame("kept.ratio" = raw.kmer.nb/nb.kmer.all, 
                                            "method" = str_replace(x, pattern = "_", replacement = ":"),
                                            "wrong.ratio" = calc_wrong_ratio(aligntab.path, fa.path),
                                            "read.type" = d))
        }
    }
}

ggplot(extens.qual) +
    geom_line(aes(x = kept.ratio, y = wrong.ratio, color = method), size = 1) +
    geom_point(aes(x = kept.ratio, y = wrong.ratio, color = method), size = 3) +
    scale_color_manual(values = c("azure4", "gold3", "royalblue", "brown")) +
    scale_x_reverse(breaks = seq(0, 1, by = 0.1),
                    labels = paste0(as.character(seq(0, 1, by = 0.1) * 100), "%")) +
    scale_y_continuous(breaks = seq(0, 0.16, by = 0.02),
                       labels = paste0(as.character(seq(0, 0.16, by = 0.02) * 100), "%")) +
    # facet_grid(read.type ~ ., scales = "free") +
    xlab("ratio of k-mers for extension") +
    ylab("ratio of mis-extension") +
    theme(text = element_text(size = 30), legend.position = "top") +
    ggsave(paste0(out.dir, "/wrong_extension_ratio_with_kept_ratio.png"), width = 12, height = 6)
 