#!/bin/bash

set -e

mkdir -p output/rankingfirst/ output/kamrat.idx

# make index, normalizating sample counts as "k-mer count per billion"
../bin/kamrat index -intab data/kmer-counts.subset4toy.tsv.gz -outdir output/kamrat.idx -klen 31 -unstrand -nfbase 1000000000
# select condition-specific k-mers
../bin/kamrat filter -idxdir output/kamrat.idx -design <(sed 's/normal/DOWN/g' data/sample-condition.toy.tsv | sed 's/tumor/UP/g') -upmin 5:5 -downmax 0:10 -outpath output/rankingfirst/top-kmers.bin
# merge top k-mers into contigs, using default intervention method (pearson:0.20), assign contig score as the minimum value among compositional k-mers, and output contig counts by calculating mean counts among compositional k-mers
../bin/kamrat merge -idxdir output/kamrat.idx -overlap 30-15 -with output/rankingfirst/top-kmers.bin:min -outpath output/rankingfirst/top-ctg-counts.tsv -withcounts mean
