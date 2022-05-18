rm(list = ls())

library(Biostrings)
library(stringr)
library(pheatmap)

workdir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"

# sample condition for heatmap annotation
sample.condi <- read.table(paste0(workdir, "/matrices/all-samples.tsv"), header = F)
names(sample.condi) <- c("sample", "condition")
sample.condi <- sample.condi[order(sample.condi$condition), ]
rownames(sample.condi) <- sample.condi$sample
    
# merged contigs
ctg.count <- read.table(paste0(workdir, "kamrat_res/all-samples/filter-merge/spec-contig-counts.tsv"),
                        header = T) # 164 contigs
ctg.annot <- read.table(paste0(workdir, "kamrat_res/all-samples/filter-merge/merged_annotation.tsv"),
                        header = T)
res <- merge(ctg.count, ctg.annot)
rownames(res) <- paste0("ctg", 1 : nrow(res), "-", res$gene_symbol)
res <- res[res$nb.merged.kmer > 15, c(-2, -3)] # 55 contigs
res$is.annot <- as.character(!is.na(res$gene_symbol))

pheatmap(res[, as.character(sample.condi$sample)],
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("white", "gray65", "black"), bias = 2)(50),
         border_color = NA,
         annotation_col = sample.condi["condition"],
         annotation_row = res["is.annot"],
         annotation_colors = list(condition = c("normal" = "royalblue", "tumor" = "gold3"),
                                  is.annot = c("TRUE" = "skyblue", "FALSE" = "tomato")),
         annotation_legend = FALSE,
         show_colnames = FALSE,
         fontsize = 20,
         cellwidth = 5,
         cellheight = 20,
         filename = paste0(workdir, "kamrat_res/all-samples/filter-merge/spec-ctg-heatmap.tiff"))
