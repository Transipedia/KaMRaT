#!/bin/bash

set -e

mkdir -p output/rankingfirst/ output/kamrat.idx

# make index, normalizating sample counts as "k-mer count per billion"
../bin/kamrat index -intab data/kmer-counts.subset4toy.tsv.gz -outdir output/kamrat.idx -klen 31 -unstrand -nfbase 1000000000
# evalate association between sample counts and sample conditions, and select top 10% of most informative contigs
../bin/kamrat rank -idxdir output/kamrat.idx -rankby ttest.padj -design data/sample-condition.toy.tsv -seltop 0.1 -outpath output/rankingfirst/top-kmers.bin
# merge top k-mers into contigs, using default intervention method (pearson:0.20), assign contig score as the minimum value among compositional k-mers, and output contig counts by calculating mean counts among compositional k-mers
../bin/kamrat merge -idxdir output/kamrat.idx -overlap 30-15 -with output/rankingfirst/top-kmers.bin:min -outpath output/rankingfirst/top-ctg-counts.tsv -withcounts mean
