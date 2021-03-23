rm(list = ls())

cmdArgs <- commandArgs(trailingOnly = T)
sig.fa.path <- cmdArgs[1]
# sig.fa.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/e_gene_level/b_rank_res/train0-lrc/signatures-lasso.fa"
count.mat.path <- cmdArgs[2]
# count.mat.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/e_gene_level/a_matrices/gene-counts.test0.tsv"
out.path <- cmdArgs[3]
# out.path <- ""

count.mat <- read.table(count.mat.path, header = T)
smp.nf <- colSums(count.mat[, -1])
smp.nf <- round(mean(smp.nf) / smp.nf, digits = 6)
for (s in names(smp.nf)) {
    count.mat[, s] <- count.mat[, s] * smp.nf[s]
}
sig.list <- readLines(sig.fa.path)
sig.list <- sig.list[seq(2, length(sig.list), by = 2)]
if (!all(sig.list %in% count.mat$feature)) {
    stop("some signatures not exist in test table")
}
write.table(count.mat[count.mat$feature %in% sig.list, ], out.path, row.names = F, col.names = T, quote = F, sep = "\t")
