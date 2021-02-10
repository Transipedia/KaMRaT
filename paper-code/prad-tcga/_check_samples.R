rm(list = ls())

library(readxl)
library(stringr)
library(magrittr)

all.smp.annot.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/data/paperHa2020-TCGAinfo.xlsx"
good.smp.annot.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/data/paperTCGA-TableS1.xls"
smp.condi.path <- paste0("/home/haoliang.xue/media/ssfa/MEMBERS/thi-ngoc-ha.nguyen/", 
                         c("Prostate_TCGA_relapse/DEkupl_run_TCGA_relapse/metadata/sample_conditions.tsv",
                           "Prostate_TCGA_LRvsHR/DEkupl_LRvsHR_limma_nonMask/Running_DEkupl/metadata/sample_conditions.tsv"))

# ===== Load Data ===== #
all.smp.annot <- read_excel(all.smp.annot.path)[, c("barcode", "filename")] %>% na.omit()
good.smp.annot <- read_excel(good.smp.annot.path, sheet = 1)[, c("SAMPLE_ID", "Reviewed_Gleason_sum", "Clinical_Gleason_sum")]
smp.condi <- lapply(smp.condi.path, FUN = function(p) read.table(file = p, header = T))

# ===== Dealing with 'NoIndex' files in smp.condi ===== #
for (i in 1 : length(smp.condi)) {
    smp.condi[[i]]$realsample <- str_replace(smp.condi[[i]]$sample, pattern = "NoIndex", replacement = "CAGATC")
}

# ===== Make all.smp.annot's filenames correspond to smp.condi ===== #
all.smp.annot$suffix <- str_extract(all.smp.annot$filename, pattern = "[0-9]+\\_[ACGT]+.tar.gz$")
all.smp.annot$numstr <- paste0("L00", str_extract(all.smp.annot$suffix, pattern = "[0-9]+"))
all.smp.annot$suffix <- str_remove(all.smp.annot$suffix, pattern = "[0-9]+_") %>%
    str_remove(pattern = ".tar.gz")
all.smp.annot$tmpfilename <- str_extract(all.smp.annot$filename, pattern = "[A-Za-z0-9\\-_]*.tar.gz$") %>%
    str_remove(pattern = "_[0-9]+\\_[ACGT]+.tar.gz$") %>%
    str_replace(pattern = "-", replacement = "_")
all.smp.annot$realfilename <- paste(all.smp.annot$tmpfilename, all.smp.annot$suffix, all.smp.annot$numstr, sep = "_")
all.smp.annot$realfilename <- paste0("X", all.smp.annot$realfilename)

# ===== Check if smp.condi names are all in all.smp.annot file ===== #
for(c in smp.condi) {
    if(!all(c$realsample %in% all.smp.annot$realfilename)) {
        stop("some sample name not in annotation")
    }
}

# ===== Make all.smp.annot's barcode correspond to good.smp.annot's SAMPLE_ID ===== #
all.smp.annot$SAMPLE_ID <- str_sub(all.smp.annot$barcode, end = 15)

# ===== Make final table by merge ===== #
final.tab <- lapply(smp.condi, 
                    FUN = function(c) merge(x = all.smp.annot[, c("barcode", "filename", "realfilename", "SAMPLE_ID")],
                                            y = c, by.x = "realfilename", by.y = "realsample") %>%
                            merge(y = good.smp.annot, by = "SAMPLE_ID"))

write.table(final.tab[[1]][, c("sample", "condition")],
            file = "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/data/sample-condition.relapse.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")
write.table(final.tab[[2]][, c("sample", "condition")],
            file = "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/PRAD_TCGA/data/sample-condition.risk.tsv",
            row.names = F, col.names = T, quote = F, sep = "\t")
