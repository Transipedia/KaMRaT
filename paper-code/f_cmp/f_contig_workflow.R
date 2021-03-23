rm(list = ls())

library()

work.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA"

getSeqRC <- function(seq)
{
    seq.rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq.rc)
}

retrieveKMers <- function(contig.list, k.len, canonical) {
    tmp.df <- NULL
    for (i in 1 : length(contig.list)) {
        ctg <- contig.list[i]
        seq.len <- nchar(ctg)
        kmer.vect <- NULL
        for (start_pos in 1 : (seq.len - k.len + 1)) {
            kmer.seq <- substr(unname(ctg), start_pos, start_pos + k.len - 1)
            if (canonical) {
                kmer.seq.rc <- getSeqRC(kmer.seq)
                if (kmer.seq.rc < kmer.seq) {
                    kmer.seq <- kmer.seq.rc
                }
            }
            kmer.vect <- c(kmer.vect, kmer.seq)
        }
        tmp.df <- rbind(tmp.df, 
                        data.frame("contig" = names(ctg), "kmer" = kmer.vect))
    }
    return(tmp.df)
}

for (d in c("risk", "relapse")) {
    for (i in 0 : 4) {
        mf.ctg <- read.table(paste0(work.dir, "/merging-first/", d, "/train", i, "/snr/top5000-contig-counts.tsv"), 
                             header = T)$contig %>% as.character()
        print(mean(nchar(mf.ctg)))
        rk.ctg <- read.table(paste0(work.dir, "/ranking-first/", d, "/train", i, "/snr/top5000-contig-counts.tsv"),
                             header = T)$contig %>% as.character()
        
        print(mean(nchar(rk.ctg)))
        mf.kmer.df <- retrieveKMers(mf.ctg, 31, T)
    }
}
