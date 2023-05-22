rm(list = ls())

library(DESeq2)

work.dir <- "../../../../working-bench/revision/benchmark-rank/"

dlist <- c("ft20000-es10-pout0", "ft20000-es1.5-pout0", "ft20000-es1.5-pout20", "ft200000-es1.5-pout20")
for (simu.set in dlist) {
    print(simu.set)
    
    cts <- read.table(paste0(work.dir, "compcodeR_res/", simu.set, "/countmat.tsv"),
                      header = TRUE, row.names = 1)
    coldata <- read.table(paste0(work.dir, "compcodeR_res/", simu.set, "/smpinfo.tsv"),
                          header = TRUE, row.names = "sample")
    coldata$condition <- factor(coldata$condition)
    all(rownames(coldata) == colnames(cts))
    
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    res.ordered <- as.data.frame(res[order(res$pvalue, -abs(res$log2FoldChange),
                                           decreasing = FALSE),])
    res.ordered$feature <- rownames(res.ordered)
    res.ordered <- res.ordered[, c("feature", "padj", "baseMean", "log2FoldChange",
                                   "lfcSE", "stat", "pvalue")]
    colnames(res.ordered)[2] <- c("deseq2.padj")
    write.table(res.ordered, row.names = FALSE, col.names = TRUE, quote = FALSE,
                file = paste0(work.dir, "kamrat_res/", simu.set, "/ranked.deseq2.padj.tsv"))
}