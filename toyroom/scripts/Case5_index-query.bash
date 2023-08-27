#!/bin/bash

set -e

BINDS='-B ../data:/sif_in -B ../output:/sif_out'
#KAMRAT='/path/to/KaMRaT.sif'
KAMRAT='../../../../software/KaMRaT.sif'

singularity exec $BINDS $KAMRAT mkdir -p /sif_out/query/ /sif_out/kamrat.idx

# make index with normalisation
singularity exec $BINDS $KAMRAT kamrat index -intab /sif_in/kmer-counts.subset4toy.tsv.gz \
					     -outdir /sif_out/kamrat.idx -klen 31 -unstrand -nfbase 1000000

# To select novel k-mers
singularity exec $BINDS $KAMRAT kamrat query -idxdir /sif_out/kamrat.idx \
					     -fasta /sif_in/sequence.toy.fa -toquery median \
					     -outpath /sif_out/query/query_counts.tsv
