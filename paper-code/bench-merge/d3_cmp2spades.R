rm(list = ls())

library(Biostrings)

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-merge/"

evaluateContigs <- function(fa.path, align.path) {
    ctg.fa <- readDNAStringSet(fa.path)
    align.res <- read.table(align.path, header = T)
    print(paste("ctg.nb:", length(ctg.fa)))
    print(paste("median.length:", median(nchar(ctg.fa))))
    align.perf.res <- align.res[align.res$qlen == align.res$qend &
                                    align.res$qstart == 1 &
                                    align.res$qlen == align.res$length &
                                    align.res$qlen == align.res$nident &
                                    align.res$pident == 100, ]
    print(paste("perfect.align.ratio:", nrow(align.perf.res) / length(ctg.fa)))
    iden.perf.res <- align.res[align.res$qlen == align.res$qend &
                                   align.res$qstart == 1 &
                                   align.res$qlen == align.res$length &
                                   align.res$qlen == align.res$nident &
                                   align.res$sstart == 1 &
                                   align.res$send == align.res$slen &
                                   align.res$slen == align.res$qlen &
                                   align.res$pident == 100, ]
    print(paste("perfect.iden.ratio:", nrow(iden.perf.res) / length(ctg.fa)))
    print(paste("perf.align and perf.identity:", sum(iden.perf.res$qseqid %in% align.perf.res$qseqid)/nrow(align.perf.res)))
}

# rnaSPAdes all reads
evaluateContigs(paste0(work.dir, "spades_res/allreads/transcripts.fasta"), 
                paste0(work.dir, "spades_res/allreads/blastn_align.tsv"))

# rnaSPAdes all k-mers
evaluateContigs(paste0(work.dir, "spades_res/allkmers/transcripts.fasta"),
                paste0(work.dir, "spades_res/allkmers/blastn_align.tsv"))

# kamrat all k-mers
evaluateContigs(paste0(work.dir, "blastn_res/100pct/ctg-seq.none.fa"),
                paste0(work.dir, "blastn_res/100pct/ctg-aligned.none.tsv"))
