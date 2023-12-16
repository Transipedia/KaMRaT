#######################################################################
# The MIT License
#
# Copyright (c) 2017, Jérôme Audoux (jerome.audoux@inserm.fr)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind,
# express or implied, including but not limited to the warranties of
# merchantability, fitness for a particular purpose and
# noninfringement. In no event shall the authors or copyright holders
# be liable for any claim, damages or other liability, whether in an
# action of contract, tort or otherwise, arising from, out of or in
# connection with the software or the use or other dealings in the
# Software.
#######################################################################

library("data.table")
library("foreach")
library("doParallel")
library("DESeq2")

args <- commandArgs(TRUE)

# Get parameters for the test
binary                    = args[1]#snakemake@input$binary
kmer_counts               = args[2]#snakemake@input$counts
sample_conditions         = args[3]#snakemake@input$sample_conditions
pvalue_threshold          = args[4]#snakemake@params$pvalue_threshold
log2fc_threshold          = args[5]#snakemake@params$log2fc_threshold
conditionA                = args[6]#snakemake@params$conditionA
conditionB                = args[7]#snakemake@params$conditionB
nb_core                   = args[8]#snakemake@threads
chunk_size                = as.numeric(args[9])#snakemake@params$chunk_size
seed                      = args[14]#snakemake@params$seed

# Get output files
output_tmp                = args[10]#snakemake@output$tmp_dir
output_diff_counts        = args[11]#snakemake@output$diff_counts
output_pvalue_all         = args[12]#snakemake@output$pvalue_all
output_log                = args[13]#snakemake@log[[1]]

# Temporary files
output_tmp_chunks         = paste(output_tmp,"/tmp_chunks/",sep="")
output_tmp_DESeq2         = paste(output_tmp,"/tmp_DESeq2/",sep="")
header_kmer_counts         = paste(output_tmp,"/header_kmer_counts.txt",sep="")
tmp_concat                = paste(output_tmp,"/tmp_concat.txt",sep="")
adj_pvalue                = paste(output_tmp,"/adj_pvalue.txt.gz",sep="")
dataDESeq2All             = paste(output_tmp,"/dataDESeq2All.txt.gz",sep="")
dataDESeq2Filtered        = paste(output_tmp,"/dataDESeq2Filtered.txt.gz",sep="")

# Create directories
dir.create(output_tmp, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_chunks, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_DESeq2, showWarnings = FALSE, recursive = TRUE)

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}
# Return the number of line in the last files of the directory
nbLineLastFile <- function(dir) {
  return(as.numeric(system(paste("zcat", paste(dir, "/$(ls ", dir, " |sort -n|grep subfile|tail -1)",sep=""), "| wc -l"), intern=TRUE)))
}
# Return the number of files in the directory
nbFiles <- function(dir) {
  return(as.numeric(system(paste("ls ", dir, "|grep subfile | wc -l"), intern=TRUE)))
}

logging("Start DESeq2_diff_methods")

# Check the chunk size
if(chunk_size > 1000000){
  logging(paste("Chunks too large for DESeq2 computations, reduce from",chunk_size,"to 1 000 000"))
  chunk_size = 1000000
}

# Set the number of cores to use
registerDoParallel(cores=nb_core)

# CLEAN THE TMP FOLDER
system(paste("rm -f ", output_tmp_chunks, "/*", sep=""))

# SAVE THE HEADER INTO A FILE
system(paste("zcat", kmer_counts, "| head -1 | cut -f2- >", header_kmer_counts))

# SHUFFLE AND SPLIT THE MAIN FILE INTO CHUNKS WITH AUTOINCREMENTED NAMES, ACCORDING TO SEED
if(seed == 'fixed'){
    cmd = paste0("zcat ", kmer_counts, " > ", output_tmp_chunks, "tmp_shuff; cat ", output_tmp_chunks, "tmp_shuff | tail -n +2 | shuf --random-source=", output_tmp_chunks, "tmp_shuff | awk -v chunk_size=", chunk_size, " -v output_tmp_chunks=", output_tmp_chunks, " 'NR % chunk_size == 1 {OFS = \"\\t\"; x = ++i\"_subfile.txt.gz\"} {OFS = \"\"; print | \"gzip >\" output_tmp_chunks x}'")
    cat(cmd, "\n")
    system(cmd)
    system(paste0("rm ", output_tmp_chunks, "/tmp_shuff"))
}else{
    system(paste("zcat", kmer_counts, "| tail -n +2 | shuf | awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
                 "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip >\" output_tmp_chunks x}'"))
}

logging("Shuffle and split done")

nb_line_last_file = nbLineLastFile(output_tmp_chunks)
nb_files = nbFiles(output_tmp_chunks)

# IF THE LAST FILE HAS LESS THAN HALF OF THE CHUNK SIZE
# CONCATENATE THE LAST 2 FILES AND THEN SPLIT
if(nb_files > 1 && nb_line_last_file < (chunk_size/2)) {

  ## CONCATENATE THE 2 FILES
  logging(paste("The last file has",nb_line_last_file,"line(s) it will be concatenated to the second last one"))

  before_last_file = paste(output_tmp_chunks, (nb_files - 1), "_subfile.txt.gz", sep="")
  last_file        = paste(output_tmp_chunks, (nb_files), "_subfile.txt.gz", sep="")

  #CONCATENATE THE LAST 2 FILES INTO A TMP FILE
  system(paste("cat", before_last_file, last_file, ">", tmp_concat, sep=" "))

  #NUMBER OF LINE OF THE TMP FILE
  nb_line_last_file = chunk_size + nb_line_last_file

  logging(paste("The last file has", nb_line_last_file, "line(s) it will be splitted in two"))

  ### DIVIDE IN TWO PARTS
  system(paste("zcat", tmp_concat, "| head -n", floor(as.integer(nb_line_last_file/2)), "| gzip >", before_last_file))
  system(paste("zcat", tmp_concat, "| tail -n", paste("+", floor(as.integer(nb_line_last_file/2 + 1)), sep=""), "| gzip >", last_file))
  system(paste("rm", tmp_concat))
}

## LOAD THE FILENAMES OF THE DIFFERENT CHUNKS
lst_files = system(paste("find",output_tmp_chunks,"-iname \"*_subfile.txt.gz\" | sort -n"), intern = TRUE)

logging("Split done")

## LOAD THE HEADER
header = as.character(unlist(read.table(file = header_kmer_counts, sep = "\t", header = FALSE)))

logging(paste("Foreach of the", length(lst_files),"files"))

## LOADING PRIOR KNOWN NORMALISATION FACTORS
colData = read.table(sample_conditions,header=T,row.names=1)

## DESeq2 ANALYSIS ON EACH CHUNKS
invisible(foreach(i=1:length(lst_files)) %dopar% {
            bigTab = read.table(lst_files[i],header=F,stringsAsFactors=F)
            #SET TAGS AS ROWNAMES
            rownames(bigTab)=bigTab[,1]
            #REMOVE THE TAG AS A COLUMN
            bigTab=bigTab[,2:ncol(bigTab)]
            names(bigTab)=header

            countData = as.matrix(bigTab)

            dds <- DESeqDataSetFromMatrix(countData,
                                          colData = colData,
                                          design = ~ condition)

            NormCount_names = colnames(bigTab)
            rm(bigTab);gc()

            #RUN DESEQ2 AND COLLECT DESeq2 results

            #REPLACE SIZE FACTORS by SIZE FACTORS COMPUTED ON
            #THE ALL DATASET
            normFactors <- matrix(colData$normalization_factor,
                                  ncol=ncol(dds),nrow=nrow(dds),
                                  dimnames=list(1:nrow(dds),1:ncol(dds)),
                                  byrow = TRUE)

            normalizationFactors(dds) <- normFactors

            #RUN DESeq2
            dds <- estimateDispersionsGeneEst(dds)

            dds <-  tryCatch(
                             estimateDispersionsMAP(estimateDispersionsFit(dds)),
                             error=function(e){
                                 cat("Error during estimateDispersionsFit, probably a '2-order-of-magnitude' message, trying suggested alterative method\n")
                                 dispersions(dds) <- mcols(dds)$dispGeneEst
                                 dds
                             })

            dds <- nbinomWaldTest(dds)
            resDESeq2 <- results(dds, pAdjustMethod = "none", contrast = c("condition",conditionB,conditionA))

            #COLLECT COUNTS
            NormCount<- as.data.frame(counts(dds, normalized=TRUE))
            names(NormCount) <- NormCount_names

            # WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
            # Kmer_ID
            # meanA
            # meanB
            # log2FC
            # NormCount

            write.table(data.frame(ID=rownames(resDESeq2),
                                   meanA=rowMeans(NormCount[,rownames(subset(colData, condition == conditionA))]),
                                   meanB=rowMeans(NormCount[,rownames(subset(colData, condition == conditionB))]),
                                   log2FC=resDESeq2$log2FoldChange,
                                   NormCount),
                        file=gzfile(paste(output_tmp_DESeq2,i,"_dataDESeq2_part_tmp.gz", sep="")),
                        sep="\t",quote=FALSE,
                        row.names = FALSE,
                        col.names = FALSE)

            # WRITE PVALUES FOR THE CURRENT CHUNK
            write.table(data.frame(ID=rownames(resDESeq2),pvalue=resDESeq2$pvalue),
                        file=gzfile(paste(output_tmp_DESeq2,i,"_pvalue_part_tmp.gz",sep="")),
                        sep="\t",quote=FALSE,
                        row.names = FALSE,
                        col.names = FALSE)

            # Remove processed chunk
            system(paste("rm",lst_files[i]))

}) #END FOREACH

system(paste("rm -rf", output_tmp_chunks))

logging("Foreach done")

#MERGE ALL CHUNKS PVALUE INTO A FILE
system(paste("find", output_tmp_DESeq2, "-name '*_pvalue_part_tmp.gz' | xargs cat >", output_pvalue_all))

logging(paste("Pvalues merged into",output_pvalue_all))

#MERGE ALL CHUNKS DESeq2 INTO A FILE
system(paste("find", output_tmp_DESeq2, "-name '*_dataDESeq2_part_tmp.gz' | xargs cat >", dataDESeq2All))

logging(paste("DESeq2 results merged into",dataDESeq2All))

# REMOVE DESeq2 CHUNKS RESULTS
system(paste("rm -rf", output_tmp_DESeq2))

#CREATE AND WRITE THE ADJUSTED PVALUE UNDER THRESHOLD WITH THEIR ID
pvalueAll         = read.table(output_pvalue_all, header=F, stringsAsFactors=F)
names(pvalueAll)  = c("tag","pvalue")
adjPvalue         = p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")

adjPvalue_dataframe = data.frame(tag=pvalueAll$tag,
                                 pvalue=adjPvalue)

write.table(adjPvalue_dataframe,
            file=gzfile(adj_pvalue),
            sep="\t",
            quote=FALSE,
            col.names = FALSE,
            row.names = FALSE)

logging("Pvalues are adjusted")

#LEFT JOIN INTO dataDESeq2All
#GET ALL THE INFORMATION (ID,MEAN_A,MEAN_B,LOG2FC,COUNTS) FOR DE KMERS
system(paste("echo \"LANG=en_EN join <(zcat ", adj_pvalue," | LANG=en_EN sort -n -k1) <(zcat ", dataDESeq2All," | LANG=en_EN sort -n -k1) | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs(\\$5) >=", log2fc_threshold, " && \\$2 <= ", pvalue_threshold, ") print \\$0}' | tr ' ' '\t' | gzip > ", dataDESeq2Filtered, "\" | bash", sep=""))
system(paste("rm", adj_pvalue, dataDESeq2All))

logging("Get counts for pvalues that passed the filter")

#CREATE THE FINAL HEADER USING ADJ_PVALUE AND DATADESeq2ALL ONES AND COMPRESS THE FILE
# CREATE THE HEADER FOR THE DESeq2 TABLE RESULT

#SAVE THE HEADER
system(paste("echo 'tag\tpvalue\tmeanA\tmeanB\tlog2FC' | paste - ", header_kmer_counts," | gzip > ", output_diff_counts))
system(paste("cat", dataDESeq2Filtered, ">>", output_diff_counts))
system(paste("rm", dataDESeq2Filtered))

logging("Analysis done")
