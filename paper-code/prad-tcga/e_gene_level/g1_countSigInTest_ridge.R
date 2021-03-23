rm(list = ls())

cmdArgs <- commandArgs(trailingOnly = T)
sig.path.prefix <- cmdArgs[1]
sig.path.prefix <- "/home/haoliang.xue/media/data/PRAD_TCGA/e_gene_level/b_rank_res/train"
test.path.prefix <- cmdArgs[2]
test.path.prefix <- "/home/haoliang.xue/media/data/PRAD_TCGA/e_gene_level/a_matrices/gene-counts.test"

for (i in 0 : 4) {
    test.path <- paste0(test.path.prefix, i, ".tsv")
    test.tab <- read.table(test.path, header = T)
    smp.nf <- colSums(test.tab[, -1])
    smp.nf <- round(mean(smp.nf) / smp.nf, digits = 6)
    for (s in names(smp.nf)) {
        test.tab[, s] <- test.tab[, s] * smp.nf[s]
    }
    for (m in c("snr")) {
        print(paste0(i, "-", m))
        sig.path <- paste0(sig.path.prefix, i, "-", m, "/signatures-ridge.fa")
        sig.list <- readLines(sig.path)
        sig.list <- sig.list[seq(2, length(sig.list), by = 2)]
        if (!all(sig.list %in% test.tab$feature)) {
            stop("some signatures not exist in test table")
        }
        out.path <- paste0(sig.path.prefix, i, "-", m, "/signatures-in-test-ridge.tsv")
        write.table(test.tab[test.tab$feature %in% sig.list, ], out.path, row.names = F, col.names = T, quote = F, sep = "\t")
    }
}
