#!/bin/bash

set -e

mkdir -p output/mergingfirst/ output/kamrat.idx

# make index, normalizating sample counts as "k-mer count per billion"
../bin/kamrat index -intab data/kmer-counts.subset4toy.tsv.gz -outdir output/kamrat.idx -klen 31 -unstrand -nfbase 1000000000
# merge k-mers into contigs, using default intervention method (pearson:0.20)
../bin/kamrat merge -idxdir output/kamrat.idx -overlap 30-15 -outpath output/mergingfirst/merged-kmers.bin
# evalate association between sample counts and sample conditions, -seltop 0.1 may be added for selecting top 10% of contigs
../bin/kamrat rank -idxdir output/kamrat.idx -rankby ttest.padj -design data/sample-condition.toy.tsv -with output/mergingfirst/merged-kmers.bin -outpath output/mergingfirst/top-ctg-counts.tsv -withcounts
