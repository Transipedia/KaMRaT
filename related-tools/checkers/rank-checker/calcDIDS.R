rm (list = ls())

res.tab <- read.table("/home/haoliang.xue/development/new-kamrat-test/rank-res/rank-diff-counts-dids.tsv",
                      header = T)
smp.info <- read.table("/home/haoliang.xue/development/new-kamrat-test/data/luad-sample-conditions.tsv", header = F)
names(smp.info) <- c("sample", "condition")
smp.normal <- as.character(smp.info$sample[smp.info$condition == "normal"])
smp.tumor <- as.character(smp.info$sample[smp.info$condition == "tumor"])

calcDIDS <- function(row) {
    norm.vect <- as.numeric(row[smp.normal])
    tumo.vect <- as.numeric(row[smp.tumor])
    norm.max <- max(norm.vect)
    tumo.max <- max(tumo.vect)
    if (tumo.max > norm.max) {
        tumo.possum <- sum(sqrt(tumo.vect[tumo.vect > norm.max] - norm.max))
        return(tumo.possum)
    } else {
        norm.possum <- sum(sqrt(norm.vect[norm.vect > tumo.max] - tumo.max))
        return(norm.possum)
    }
}

dids.R <- apply(res.tab, MARGIN = 1, function(r) calcDIDS(r))
print(max(abs(dids.R - res.tab$DIDS.score) / dids.R))
