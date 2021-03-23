rm(list = ls())

library(rjson)
library(magrittr)
library(ROCR)

dir.prefix <- "/home/haoliang.xue/media/data/PRAD_TCGA_relapse/e_imoka/c_testing_res/train"

for (i in 0 : 4) {
    smp.info <- read.table(paste0(dir.prefix, i, "/matrix.test.tsv"), header = T, nrows = 1) %>%
        t() %>% 
        as.data.frame()
    
    json.content <- fromJSON(file = paste0(dir.prefix, i, "/predictions.tsv"))
    pred.prob <- json.content$probabilities %>% 
        do.call(what = rbind) %>% 
        as.data.frame()
    colnames(pred.prob) <- c("prob0", "prob1")
    pred.prob$sample <- json.content["samples"] %>% unlist()
    
    if (!all(rownames(smp.info) == pred.prob$sample)) {
        stop("sample info not coherent with pred prob !")
    }
    pred4auc <- prediction(predictions = pred.prob[, 2], 
                           labels = ifelse(as.character(smp.info$group) == "Yes", yes = 1, no = 0))
    auc.pred <- performance(pred4auc, measure = "auc")@y.values[[1]]
    print(paste0("AUC of train", i, ": ", round(auc.pred, 3)))
}