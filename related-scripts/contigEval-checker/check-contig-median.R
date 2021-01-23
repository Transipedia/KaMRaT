rm(list = ls())

cmdArgs <- commandArgs(trailingOnly = T)
kamrat.tab.path <- cmdArgs[1]
# kamrat.tab.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/count-in-tests/contig-counts.test0.tsv"
r.tab.path <- cmdArgs[2]
# r.tab.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-urine/cv.results/b_deseq2_contigs/count-in-tests/contig-counts.test0.byR.tsv"

kamrat.tab <- read.table(kamrat.tab.path, header = T, row.names = 1)
r.tab <- read.table(r.tab.path, header = T, row.names = 1)

diff.tab <- abs(kamrat.tab - r.tab)
for (i in 1 : nrow(diff.tab)) {
    for (j in 1 : ncol(diff.tab)) {
        if (r.tab[i, j] != kamrat.tab[i, j]) {
            diff.tab[i, j] <- diff.tab[i, j] / r.tab[i, j]
        }
    }
}
max(diff.tab)
