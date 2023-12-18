#!/bin/bash

set -e

BINDS='-B ../data:/sif_in -B ../output:/sif_out'
#KAMRAT='/path/to/KaMRaT.sif'
KAMRAT='../../../../software/KaMRaT.sif'

singularity exec $BINDS $KAMRAT mkdir -p /sif_out/rank-merge/ /sif_out/kamrat.idx

# make index with normalisation
singularity exec $BINDS $KAMRAT kamrat index -intab /sif_in/kmer-counts.subset4toy.tsv.gz \
       					     -outdir /sif_out/kamrat.idx -klen 31 -unstrand -nfbase 1000000

# To evalate association between sample counts and sample states by t-test,
# and select top 10% of most informative contigs
singularity exec $BINDS $KAMRAT kamrat rank -idxdir /sif_out/kamrat.idx -rankby ttest.padj \
					    -design /sif_in/sample-states.toy.tsv -seltop 0.1 \
					    -outpath /sif_out/rank-merge/top-kmers.bin

# To merge top k-mers into contigs, using default intervention method (pearson:0.20),
# assign contig score as the minimum value among compositional k-mers, 
# and output contig counts by calculating mean counts among compositional k-mers
singularity exec $BINDS $KAMRAT kamrat merge -idxdir /sif_out/kamrat.idx -overlap 30-15 \
					     -with /sif_out/rank-merge/top-kmers.bin:min \
					     -outpath /sif_out/rank-merge/top-ctg-counts.tsv -withcounts mean
