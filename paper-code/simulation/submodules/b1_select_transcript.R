rm(list = ls())

library(Biostrings)
library(magrittr)
library(stringr)

cmdArgs <- commandArgs(trailingOnly = TRUE)
fasta.path <- cmdArgs[1]
out.path <- cmdArgs[2]

# fasta.path <- "/home/haoliang.xue/media/data/kamrat/paper/a_ensembl_data/gencode.v34.transcripts.fa"
# out.path <- "/home/haoliang.xue/media/data/kamrat/paper/b_reference_gen/a_selected_transcripts.fa"

print(paste("fasta path:", fasta.path))
print(paste("output path:", out.path))

set.seed(91400)

tx.all <- readDNAStringSet(fasta.path)
tx.all <- tx.all[nchar(tx.all) > 100] # will simulate reads with length as 100

gene.all <- names(tx.all) %>%
    str_extract(pattern = "\\|ENSG.*?\\|") %>%
    unique() # note: gene ENSEMBL ID may contain "_PAR_Y"

gene.sel <- sample(gene.all, size = 1000, replace = F) # select randomly 1000 genes
tx.sel <- tx.all[str_extract(names(tx.all), pattern = "\\|ENSG.*?\\|") %in% gene.sel] # keep all transcripts related to the selected 1000 genes

writeXStringSet(tx.sel, filepath = out.path)
