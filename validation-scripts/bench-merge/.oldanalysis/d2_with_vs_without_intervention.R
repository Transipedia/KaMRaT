rm(list = ls())

library(Biostrings)
library(tidyr)
library(ggplot2)


work.dir <- "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
# work.dir <- "../../../../Ariticles/KaMRaT/RevisedAnalysis/bench-merge/"

min_len <- 61

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

pdf(paste0(work.dir, "results/correctness_assembly_extension-minlen", min_len, ".pdf"),
    width=16, height=10)
for (rec_abd in list(c(1, 1), c(2, 1),  c(2, 3), c(2, 5), c(5, 1), c(5, 3))) {
    stats.res <- NULL
    for (dpt in c(0.1, 0.2, 0.5, 1, 2, 5, 10)) {
        # KaMRaT none
        cat(dpt, "KaMRaT none", "\n")
        ctg.path <- paste0(work.dir, "kamrat_res-",
			   rec_abd[1], "-", rec_abd[2], "/depth_", dpt, "/ctg-seq.none.fa")
        align.path <- paste0(work.dir, "kamrat_res-",
			     rec_abd[1], "-", rec_abd[2], "/depth_", dpt, "/ctg-aligned.none.tsv")
        ctg.fa <- readDNAStringSet(ctg.path)
        ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
        align.res <- read.table(align.path, header = TRUE, row.names = 1)
        align.res <- align.res[align.res$qlen > min_len, ]
        if (!all(nchar(ctg.fa) == align.res$qlen[rownames(ctg.fa)])) {
            stop(paste(ctg.path, align.path, sep = "\n"))
        }
        stats.res <- rbind(stats.res,
                           data.frame("depth" = dpt,
                                      "mode" = "KaMRaT none",
            			  "nb.ctg" = length(ctg.fa),
            			  "ctg.median.len" = median(nchar(ctg.fa)),
                                      "perf.align" = evaluate_local(ctg.fa, align.res),
            			  "identical" = evaluate_global(ctg.fa, align.res)))
        # KaMRaT with intervention
        for (mode in c("pearson", "spearman", "mac")) {
            cat(dpt, "KaMRaT", paste(mode, 0.2, sep = ":"), "\n")
            ctg.path <- paste0(work.dir, "kamrat_res-",
			       rec_abd[1], "-", rec_abd[2], "/depth_", dpt, "/ctg-seq.", mode, "_0.2", ".fa")
            align.path <- paste0(work.dir, "kamrat_res-",
				 rec_abd[1], "-", rec_abd[2], "/depth_", dpt, "/ctg-aligned.", mode, "_0.2", ".tsv")
            ctg.fa <- readDNAStringSet(ctg.path)
            ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
            align.res <- read.table(align.path, header = TRUE, row.names = 1)
            align.res <- align.res[align.res$qlen > min_len, ]
            if (!all(nchar(ctg.fa) == align.res$qlen[rownames(ctg.fa)])) {
                stop(paste(ctg.path, align.path, sep = "\n"))
            }
            stats.res <- rbind(stats.res,
                               data.frame("depth" = dpt,
                                          "mode" = paste0("KaMRaT ", mode, ":0.2"),
            		   	          "nb.ctg" = length(ctg.fa),
            			          "ctg.median.len" = median(nchar(ctg.fa)),
                                          "perf.align" = evaluate_local(ctg.fa, align.res),
            			          "identical" = evaluate_global(ctg.fa, align.res)))
        }
        # rnaSPAdes
        for (mode in c("allreads", paste0("allkmers-", rec_abd[1], "-", rec_abd[2]))) {
            cat(dpt, "SPAdes", mode, "\n")
            ctg.path <- paste0(work.dir, "spades_res/depth_", dpt, "/", mode, "/transcripts.fasta")
            align.path <- paste0(work.dir, "spades_res/depth_", dpt, "/", mode, "/blastn_align.tsv")
            ctg.fa <- readDNAStringSet(ctg.path)
      	    ctg.fa <- ctg.fa[width(ctg.fa) > min_len]
            align.res <- read.table(align.path, header = TRUE, row.names = 1)
            align.res <- align.res[align.res$qlen > min_len, ]
            if (!all(nchar(ctg.fa) == align.res$qlen[rownames(ctg.fa)])) {
                stop(paste(ctg.path, align.path, sep = "\n"))
            }
            stats.res <- rbind(stats.res,
                               data.frame("depth" = dpt,
                                          "mode" = paste0("rnaSPAdes ", strsplit(mode, split = "-")[[1]][1]),
            		   	          "nb.ctg" = length(ctg.fa),
            			          "ctg.median.len" = median(nchar(ctg.fa)),
                                          "perf.align" = evaluate_local(ctg.fa, align.res),
            			          "identical" = evaluate_global(ctg.fa, align.res)))
        }
    }
    write.csv(stats.res, paste0(work.dir, "results/cmp_assembly_extension-",
   			        rec_abd[1], "-", rec_abd[2], "-minlen", min_len, ".csv"),
	      quote = FALSE)
    plt <- ggplot(data = stats.res) +
        geom_line(aes(x = depth, y = perf.align, color = mode), linewidth = 1) +
        geom_point(aes(x = depth, y = perf.align, color = mode)) +
        scale_x_log10(breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10)) +
        scale_color_manual(values = c("KaMRaT none" = "#e66101",
	   			      "KaMRaT mac:0.2" = "#fdb863", 
				      "KaMRaT pearson:0.2" = "#b2abd2",
				      "KaMRaT spearman:0.2" = "#5e3c99",
				      "rnaSPAdes allkmers" = "#808080",
				      "rnaSPAdes allreads" = "#000000")) +
        ylab("%perfect alignment") +
	ggtitle(paste0("min_rec = ", rec_abd[1], ", min_abd = ", rec_abd[2])) +
        theme_light() +
        theme(text = element_text(size = 30, family = "sans"))
    plot(plt)
}
dev.off()
