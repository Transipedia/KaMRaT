rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)
library(ggrepel)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge-new/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

evaluate_local <- function(ctg.fa, align.res) {
    perfect.align <- (align.res$qlen == align.res$qend &
                      align.res$qstart == 1 &
                      align.res$qlen == align.res$length &
                      align.res$qlen == align.res$nident &
                      align.res$pident == 100)
    return(round(sum(perfect.align) / length(ctg.fa) * 100, 2))
}

evaluate_global <- function(ctg.fa, align.res) {
    ident <- (align.res$qlen == align.res$qend &
              align.res$qstart == 1 &
              align.res$qlen == align.res$length &
              align.res$qlen == align.res$nident &
              align.res$sstart == 1 &
              align.res$send == align.res$slen &
              align.res$slen == align.res$qlen &
              align.res$pident == 100)
    return(round(sum(ident) / length(ctg.fa) * 100, 2))
}

stats.res <- NULL
for (dpt in c(0.1)) {
    # KaMRaT none
    cat(dpt, "KaMRaT none", ":")
    ctg.path <- paste0(work.dir, "kamrat_res_err-free_3-1/depth_", dpt, "/ctg-seq.none.fa")
    align.path <- paste0(work.dir, "kamrat_res_err-free_3-1/depth_", dpt, "/ctg-aligned.none.tsv")
    ctg.fa <- readDNAStringSet(ctg.path)
    align.res <- read.table(align.path, header = TRUE, row.names = 1)
    cat("\t", length(ctg.fa), nrow(align.res), "\n")
    stats.res <- rbind(stats.res,
                       data.frame("depth" = dpt,
                                  "mode" = "KaMRaT none",
        			  "nb.ctg" = length(ctg.fa),
        			  "ctg.median.len" = median(nchar(ctg.fa)),
                                  "perf.align" = evaluate_local(ctg.fa, align.res),
        			  "identical" = evaluate_global(ctg.fa, align.res)))
    # rnaSPAdes
#    for (mode in c("allreads", "allkmers-3-1")) {
#        cat(dpt, "SPAdes", mode, ":")
#        ctg.path <- paste0(work.dir, "spades_res/err-free/depth_", dpt, "/", mode, "/transcripts.fasta")
#        align.path <- paste0(work.dir, "spades_res/err-free/depth_", dpt, "/", mode, "/blastn_align.tsv")
#        ctg.fa <- readDNAStringSet(ctg.path)
#        align.res <- read.table(align.path, header = TRUE, row.names = 1)
#	cat("\t", length(ctg.fa), nrow(align.res), "\n")
#        stats.res <- rbind(stats.res,
#                           data.frame("depth" = dpt,
#                                      "mode" = paste0("rnaSPAdes ", strsplit(mode, split = "-")[[1]][1]),    
#				      "nb.ctg" = length(ctg.fa),
#				      "ctg.median.len" = median(nchar(ctg.fa)),
#                                      "perf.align" = evaluate_local(ctg.fa, align.res),
#				      "identical" = evaluate_global(ctg.fa, align.res)))
#    }
}
write.csv(stats.res, paste0(work.dir, "results/1_cmp_assembly_extension.csv"),
          quote = FALSE)

plt <- ggplot(data = stats.res, aes(x = depth, y = perf.align, color = mode)) +
    geom_line(linewidth = 1) +
    geom_point() +
    # geom_text_repel(aes(label = ctg.median.len)) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.5, 1, 5, 10, 20, 30, 40, 50)) +
    scale_color_manual(values = c("KaMRaT none" = "#e66101",
				  "rnaSPAdes allkmers" = "#808080",
				  "rnaSPAdes allreads" = "#000000")) +
    ylab("%perfect alignment") +
    theme_light() +
    theme(text = element_text(size = 15, family = "sans"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(plt, filename = paste0(work.dir, "results/1_cmp_assembly_extension.pdf"),
       width = 9, height = 3)
