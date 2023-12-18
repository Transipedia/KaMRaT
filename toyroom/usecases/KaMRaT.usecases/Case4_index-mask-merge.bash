#!/bin/bash

set -e

BINDS='-B ../data:/sif_in -B ../output:/sif_out'
#KAMRAT='/path/to/KaMRaT.sif'
KAMRAT='../../../../software/KaMRaT.sif'

singularity exec $BINDS $KAMRAT mkdir -p /sif_out/mask-merge/ /sif_out/kamrat.idx

# make index with normalisation
singularity exec $BINDS $KAMRAT kamrat index -intab /sif_in/kmer-counts.subset4toy.tsv.gz \
					     -outdir /sif_out/kamrat.idx -klen 31 -unstrand -nfbase 1000000

# To select novel k-mers
singularity exec $BINDS $KAMRAT kamrat mask -idxdir /sif_out/kamrat.idx \
					    -fasta /sif_in/sequence.toy.fa  \
					    -outpath /sif_out/mask-merge/masked-kmers.bin

# To merge top k-mers into contigs, using default intervention method (pearson:0.20),
# assign contig score as the minimum value among compositional k-mers,
# and output contig counts by calculating mean counts among compositional k-mers
singularity exec $BINDS $KAMRAT kamrat merge -idxdir /sif_out/kamrat.idx -overlap 30-15 \
					     -with /sif_out/mask-merge/masked-kmers.bin:min \
					     -outpath /sif_out/mask-merge/masked-ctg-counts.tsv -withcounts mean
