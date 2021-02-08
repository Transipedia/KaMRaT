rm(list = ls())

library(magrittr)
library(Biostrings)
library(stringr)
library(ggplot2)
library(tidyr)

cmdArgs <- commandArgs(trailingOnly = T)
kamrat.res.dir <- cmdArgs[1]
spades.res.dir <- cmdArgs[2]
out.dir <- cmdArgs[3]

# kamrat.res.dir <- "/home/haoliang.xue/media/data/kamrat/paper/f_kamrat_res"
# spades.res.dir <- "/home/haoliang.xue/media/data/kamrat/paper/d_SPAdes_res"
# out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/g_KaMRaT_merge_eval"

get_N50 <- function(contig.list) {
    len.list <- nchar(contig.list)
    len.list <- len.list[order(len.list, decreasing = T)]
    len.sum <- sum(len.list)
    sum.x <- 0
    for (i in 1 : length(len.list)) {
        sum.x <- sum.x + len.list[i]
        if (sum.x >= len.sum / 2) {
            break
        }
    }
    return(len.list[i])
}

res.summary <- NULL

# Error-free
spades.res.dir1 <- paste0(spades.res.dir, "/a_errfree")
spades.n50 <- NULL
for (s in dir(spades.res.dir1)) {
    fname <- paste0(spades.res.dir1, "/", s, "/transcripts.fasta")
    spades.n50 <- c(spades.n50, get_N50(readDNAStringSet(fname)))
}
res.summary <- rbind(res.summary,
                     data.frame("software" = "SPAdes",
                                "N50" = paste0(min(spades.n50), "-", max(spades.n50)),
                                "dataset" = "error-free"))
rm(fname, spades.res.dir1, spades.n50, s)
kamrat.res.prefix <- paste0(kamrat.res.dir, "/a_errfree_randsel/merged-contigs-100pct-1-1.")
for (m in c("mac_0.32", "none", "pearson_0.21", "spearman_0.23")) {
    fname <- paste0(kamrat.res.prefix, m, ".fa")
    res.summary <- rbind(res.summary,
                         data.frame("software" = paste0("KaMRaT-", str_extract(m, pattern = "[a-z]*")), 
                                    "N50" = as.character(get_N50(readDNAStringSet(fname))),
                                    "dataset" = "error-free"))
}
rm(fname, kamrat.res.prefix, m)

# With-error
spades.res.dir1 <- paste0(spades.res.dir, "/b_witherr")
spades.n50 <- NULL
for (s in dir(spades.res.dir1)) {
    fname <- paste0(spades.res.dir1, "/", s, "/transcripts.fasta")
    spades.n50 <- c(spades.n50, get_N50(readDNAStringSet(fname)))
}
res.summary <- rbind(res.summary,
                     data.frame("software" = "SPAdes",
                                "N50" = paste0(min(spades.n50), "-", max(spades.n50)),
                                "dataset" = "with-error"))
rm(fname, spades.res.dir1, spades.n50, s)
kamrat.res.prefix <- paste0(kamrat.res.dir, "/c_witherr_filter/merged-contigs-100pct")
for (f in c("-1-1", "-1-2", "-1-3", "-1-4", "-1-5", "-3-5")) {
    for (m in c("mac_0.32", "none", "pearson_0.21", "spearman_0.23")) {
        print(paste0(f, "-", m))
        fname <- paste0(kamrat.res.prefix, f, ".", m, ".fa")
        res.summary <- rbind(res.summary,
                             data.frame("software" = paste0("KaMRaT-", str_extract(m, pattern = "[a-z]*")), 
                                        "N50" = as.character(get_N50(readDNAStringSet(fname))),
                                        "dataset" = paste0("with-error", f)))
    }
}
rm(fname, kamrat.res.prefix, m, f)

res.summary.wide <- pivot_wider(res.summary, id_cols = "software", names_from = "dataset", values_from = "N50")
write.table(res.summary.wide, paste0(out.dir, "/N50_comparison.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
