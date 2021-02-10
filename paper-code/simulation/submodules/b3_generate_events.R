rm(list = ls())

library(stringr)
library(magrittr)
library(httr)
library(jsonlite)
library(xml2)

cmdArgs <- commandArgs(trailingOnly = T)
gvf.path <- cmdArgs[1]
out.path <- cmdArgs[2]

# gvf.path <- "/home/haoliang.xue/media/data/kamrat/paper/.old.res/b_reference_gen/b_1000GENOMES-phase_3.sel.gvf"
# out.path <- "/home/haoliang.xue/media/data/kamrat/paper/.old.res/b_reference_gen/c_variation_events.fa"

gvf.tab <- read.table(gvf.path, header = F, sep = "\t")[, c(1, 4, 5, 7, 3, 9)]
names(gvf.tab) <- c("chr", "start", "end", "strand", "event_type", "info")
gvf.tab$strand <- ifelse(gvf.tab$strand == '+', 1, -1)
gvf.tab$var.id <- str_extract(gvf.tab$info, pattern = "Dbxref=[A-Za-z0-9_]+:[a-z0-9]+") %>% 
    str_extract(pattern = "[a-z0-9]+$")
gvf.tab$ref.seq <- str_extract(gvf.tab$info, pattern = "Reference_seq=[ACGT-]+") %>% str_extract(pattern = "[ACGT-]+")
gvf.tab$var.seq.list <- str_extract(gvf.tab$info, pattern = "Variant_seq=[ACGT,-]+") %>% str_extract_all(pattern = "[ACGT-]+")

server <- "http://rest.ensembl.org"
ext.list <- paste0("/sequence/region/human/", 
                   gvf.tab$chr, ":", gvf.tab$start - 99, "..", gvf.tab$end + 99, ":", gvf.tab$strand)

gvf.tab$seq.with.flank <- lapply(ext.list, 
                                 function(ext) {
                                     rx <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
                                     stop_for_status(rx)
                                     return(content(rx))}) %>% as.character()

event.fa <- NULL
for (i in 1 : nrow(gvf.tab)) {
    # Sequence in Reference
    ref.seq <- gvf.tab$ref.seq[i]
    if ((ref.seq == "-" && gvf.tab$end[i] != gvf.tab$start[i]) ||
        (ref.seq != "-" && sum(strsplit(ref.seq, NULL)[[1]] %in% c("A", "C", "G", "T")) != gvf.tab$end[i] - gvf.tab$start[i] + 1)) {
        stop("Length between variant zone and reference sequence not identical")
    }
    ref.seq.with.flank <- gvf.tab$seq.with.flank[i]
    if (ref.seq != "-" && str_sub(ref.seq.with.flank, start = 100, end = -100) != ref.seq) {
        stop(paste("Flanking sequence not coherent with event sequence:", 
                   str_sub(ref.seq.with.flank, start = 100, end = -100), 
                   "vs", 
                   ref.seq))
    }
    left.flank <- str_sub(ref.seq.with.flank, end = 99)
    right.flank <- ifelse(ref.seq == "-", str_sub(ref.seq.with.flank, start = -100), str_sub(ref.seq.with.flank, start = -99))
    # Sequence as Variant
    for (j in 1 : length(gvf.tab$var.seq.list[i][[1]])) {
        var.seq <- gvf.tab$var.seq.list[i][[1]][j]
        event.header <- paste0(">", gvf.tab$var.id[i], "|", gvf.tab$event_type[i], "|var", j)
        if (var.seq == "-") {
            var.seq.with.flank <- paste0(left.flank, right.flank)
        } else {
            var.seq.with.flank <- paste0(left.flank, var.seq, right.flank)
        }
        event.fa <- c(event.fa, event.header, var.seq.with.flank)
    }
}
print(any(table(event.fa[seq(2, length(event.fa), 2)]) > 1))

out.fa <- file(out.path, "w")
writeLines(event.fa, out.fa)
close(out.fa)
