library(magrittr)

getSeqRC <- function(seq)
{
    seq.rc <- paste(rev(unlist(strsplit(chartr("ACGT", "TGCA", seq), NULL))), collapse = "")
    return(seq.rc)
}

extractKMerVect <- function(seq, k.len, is_stranded) {
    kmer.vect <- c()
    for (i in 1 : (nchar(seq) - k.len + 1)) {
        kmer <- substr(seq, start = i, stop = i + k.len - 1)
        if (!is_stranded) {
            kmer.rc <- getSeqRC(kmer)
            if (kmer.rc < kmer) {
                kmer <- kmer.rc
            }
        }
        kmer.vect <- c(kmer.vect, kmer)
    }
    return(kmer.vect)
}

cmpAcrossContigLists <- function(lcl, is_stranded, k.len, show.unique) {
    ctg.cmp.df <- lapply(lcl[[1]], 
                         FUN = function(s) return(data.frame("kmer" = extractKMerVect(seq = as.character(s),
                                                                                      k.len = k.len, 
                                                                                      is_stranded = is_stranded),
                                                             "contig" = as.character(s)))) %>%
    
        do.call(what = rbind)
    names(ctg.cmp.df) <- c("kmer", names(lcl)[1])
    for (i in 2 : length(lcl)) {
        ctg.i.df <- lapply(lcl[[i]], 
                           FUN = function(s) return(data.frame("kmer" = extractKMerVect(seq = as.character(s),
                                                                                        k.len = k.len, 
                                                                                        is_stranded = is_stranded),
                                                               "contig" = as.character(s)))) %>%
            do.call(what = rbind)
        names(ctg.i.df) <- c("kmer", names(lcl)[i])
        ctg.cmp.df <- merge(x = ctg.cmp.df, y = ctg.i.df, by = "kmer", all = show.unique)
    }
    if (nrow(ctg.cmp.df) == 0) {
        return(ctg.cmp.df)
    }
    for (i in 2 : (length(lcl) + 1)) {
        ctg.cmp.df[, i] <- factor(ctg.cmp.df[, i], levels = c(levels(ctg.cmp.df[, i]), "-"))
        ctg.cmp.df[is.na(ctg.cmp.df[, i]), i] <- "-"
    }
    ctg.cmp <- aggregate(formula = kmer ~ ., data = ctg.cmp.df, FUN = length)
    names(ctg.cmp) <- c(names(lcl), "nb.kmer")
    ctg.cmp$nb.shared.contig <- apply(ctg.cmp, MARGIN = 1, function(r) return(sum(r != "-") - 1))
    ctg.cmp <- ctg.cmp[order(ctg.cmp$nb.shared.contig, ctg.cmp$nb.kmer, decreasing = T), ]
    return(ctg.cmp)
}
