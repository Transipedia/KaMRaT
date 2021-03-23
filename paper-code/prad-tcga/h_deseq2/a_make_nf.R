rm(list = ls())

dkpl.path <- "/home/haoliang.xue/media/ssfa/MEMBERS/thi-ngoc-ha.nguyen/Prostate_TCGA_LRvsHR/DEkupl_LRvsHR_limma_nonMask/Running_DEkupl/metadata/sample_conditions_full.tsv"
proj.dir <- "/home/haoliang.xue/media/ssfa/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga"

dkpl.smp.nf <- read.table(dkpl.path, header = T)

for (i in 0 : 4) {
    smp.condi <- read.table(paste0(proj.dir, "/b_splitCV/sampleshuf.train", i, ".tsv"), header = T)
    smp.condi.full <- merge(smp.condi, dkpl.smp.nf, by = c("sample", "condition"), sort = F)
    print(all(smp.condi.full$sample == smp.condi$sample))
    write.table(smp.condi.full, paste0(proj.dir, "/h_mergef-deseq2/data/smp_condi_full.train", i, ".tsv"), 
                row.names = F, col.names = T, quote = F, sep = "\t")
}
