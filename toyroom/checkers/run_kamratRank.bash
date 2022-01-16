#!/bin/bash

set -e

cd ../../

mkdir -p toyroom/output/index/ toyroom/output/kamrat-rank/
bin/kamrat index -intab toyroom/data/kmer-counts.subset4toy.tsv.gz -outdir toyroom/output/index -klen 31 -unstrand -nfbase 10000000

for m in `echo ttest.padj ttest.pi snr dids lr bayes`
do
	bin/kamrat rank -idxdir toyroom/output/index/ -rankby $m -design toyroom/data/sample-condition.toy.tsv -outpath toyroom/output/kamrat-rank/kamrat-rank-${m/:/}.txt -withcounts
done

for m in `echo pearson spearman`
do
	bin/kamrat rank -idxdir toyroom/output/index/ -rankby $m -design toyroom/data/sample-condition.toy2.tsv -outpath toyroom/output/kamrat-rank/kamrat-rank-$m.txt -withcounts
done

for m in `echo sd rsd1 rsd2 entropy`
do
	bin/kamrat rank -idxdir toyroom/output/index/ -rankby $m -outpath toyroom/output/kamrat-rank/kamrat-rank-$m.txt -withcounts
done

rm -rf toyroom/output/index/
