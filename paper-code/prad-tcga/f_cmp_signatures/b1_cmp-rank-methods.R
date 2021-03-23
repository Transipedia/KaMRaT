rm(list = ls())

library(Biostrings)
library(stringr)
library(UpSetR)

cmdArgs <- commandArgs(trailingOnly = T)
out.path <- cmdArgs[1]
# out.path <- "/home/haoliang.xue/media/data/PRAD_TCGA/f_signature_cmp/b_rank-methods/upset.png"
sig.path.list <- NULL
for (i in 2 : length(cmdArgs)) {
    sig.path.list <- c(sig.path.list, cmdArgs[i])
    print(paste0("fa", i - 2, ": ", cmdArgs[i]))
}
# sig.path.list <- paste0("/home/haoliang.xue/media/data/PRAD_TCGA/d_signatures/d4_rank-methods/model.train0-",
#                         c("lrc", "nbc", "svm", "snr", "ttest"), "/signatures-lasso.fa")

sig.list <- lapply(sig.path.list, function(sf) readDNAStringSet(sf) %>% as.character() %>% unname())
names(sig.list) <- str_extract(sig.path.list, pattern = "snr|ttest|lrc|nbc|svm")

png(out.path, width = 800, height = 500)
upset(data = fromList(sig.list), order.by = "freq", nintersects = NA, text.scale = 3)
dev.off()
