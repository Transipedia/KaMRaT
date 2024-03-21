#!/bin/bash

set -e

BINDS='-B ../data:/sif_in -B ../output:/sif_out'
#KAMRAT='/path/to/KaMRaT.sif'
KAMRAT='../../../../software/KaMRaT.sif'

singularity exec $BINDS $KAMRAT mkdir -p /sif_out/merge-rank/ /sif_out/kamrat.idx/

# make index with normalisation
singularity exec $BINDS $KAMRAT kamrat index -intab /sif_in/kmer-counts.subset4toy.tsv.gz \
					     -outdir /sif_out/kamrat.idx -klen 31 -unstrand -nfbase 1000000

# merge k-mers into contigs, using default intervention method (pearson:0.20)
singularity exec $BINDS $KAMRAT kamrat merge -idxdir /sif_out/kamrat.idx -overlap 30-15 \
					     -outpath /sif_out/merge-rank/merged-kmers.bin

# Case 1-1: To evalate association between sample counts and sample states using t-test
singularity exec $BINDS $KAMRAT kamrat rank -idxdir /sif_out/kamrat.idx -scoreby ttest.padj \
					    -design /sif_in/sample-states.toy.tsv \
					    -with /sif_out/merge-rank/merged-kmers.bin \
					    -outpath /sif_out/merge-rank/top-ctg-counts.ttest.padj.tsv -withcounts

# Case 1-2: To evaluate Spearman correlation between sample counts and sample values
singularity exec $BINDS $KAMRAT kamrat rank -idxdir /sif_out/kamrat.idx -scoreby spearman \
					    -design /sif_in/sample-values.toy.tsv \
					    -with /sif_out/merge-rank/merged-kmers.bin \
					    -outpath /sif_out/merge-rank/top-ctg-counts.spearman.tsv -withcounts

# Case 1-3: To evaluate sample counts' dispersion by standard deviation
singularity exec $BINDS $KAMRAT kamrat rank -idxdir /sif_out/kamrat.idx -scoreby sd \
					    -with /sif_out/merge-rank/merged-kmers.bin \
					    -outpath /sif_out/merge-rank/top-ctg-counts.sd.tsv -withcounts
