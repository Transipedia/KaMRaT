rm(list = ls())

cmdArgs <- commandArgs(trailingOnly = T)
list.byR.path <- cmdArgs[1]
list.kamrat.path <- cmdArgs[2]

# list.byR.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-test/norm/tmp/masked-counts.rankSNR.byR.list"
# list.kamrat.path <- "/home/haoliang.xue/media/data/kamrat/kamrat-test/norm/tmp/masked-counts.rankSNR.list"

list.byR <- read.table(list.byR.path, header = T)
names(list.byR) <- c("tag", "score.byR")
list.kamrat <- read.table(list.kamrat.path, header = T)
names(list.kamrat) <- c("tag", "score.kamrat")

list.cmp <- merge(list.byR, list.kamrat, by = "tag", all = T)

print(paste("Max of absolute difference:", max(abs(list.cmp$score.byR - list.cmp$score.kamrat))))
