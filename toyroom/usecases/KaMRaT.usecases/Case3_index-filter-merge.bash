#!/bin/bash

set -e

BINDS='-B ../data:/sif_in -B ../output:/sif_out'
#KAMRAT='/path/to/KaMRaT.sif'
KAMRAT='../../../../software/KaMRaT.sif'

singularity exec $BINDS $KAMRAT mkdir -p /sif_out/filter-merge/ /sif_out/kamrat.idx

# make index with normalisation
singularity exec $BINDS $KAMRAT kamrat index -intab /sif_in/kmer-counts.subset4toy.tsv.gz \
					     -outdir /sif_out/kamrat.idx -klen 31 -unstrand -nfbase 1000000

# To select condition-specific k-mers
singularity exec $BINDS $KAMRAT kamrat filter -idxdir /sif_out/kamrat.idx \
					      -design /sif_in/sample-indications.toy.tsv \
					      -upmin 3:5 -downmax 0:10 -outpath /sif_out/filter-merge/top-kmers.bin

# To merge top k-mers into contigs, using default intervention method (pearson:0.20),
# assign contig score as the minimum value among compositional k-mers,
# and output contig counts by calculating mean counts among compositional k-mers
singularity exec $BINDS $KAMRAT kamrat merge -idxdir /sif_out/kamrat.idx -overlap 30-15 \
					     -with /sif_out/filter-merge/filtered-kmers.bin:min \
					     -outpath /sif_out/filter-merge/filtered-ctg-counts.tsv -withcounts mean
