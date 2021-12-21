rm(list = ls())

library(Biostrings)
library(stringr)
library(pheatmap)

workdir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"

# merged contigs
ctg.annot <- read.table(paste0(workdir, "kamrat_res/all-samples/filter-merge/merged_annotation.tsv"),
                        header = T)
star.missed <- as.character(ctg.annot$contig[is.na(ctg.annot$gene_symbol)])

blastn.align <- read.table(paste0(workdir, "kamrat_res/all-samples/filter-merge/spec-contig-align.tsv"),
                           header = T)
ctg.fa <- readDNAStringSet(paste0(workdir, "kamrat_res/all-samples/filter-merge/spec-contig-seq.fa"))
ctg.seq <- data.frame("qseqid" = names(ctg.fa), "contig" = as.character(ctg.fa))

blastn.align <- merge(ctg.seq, blastn.align)
blastn.align <- blastn.align[blastn.align$contig %in% star.missed, ]
blastn.align$palign <- blastn.align$nident / blastn.align$qlen

sum(blastn.align$palign <= 0.5)

write.table(blastn.align[order(blastn.align$palign), ],
            file = paste0(workdir, "kamrat_res/all-samples/filter-merge/cmp-star-blastn.tsv"),
            row.names = F, col.names = T, quote = F, sep = "\t")
