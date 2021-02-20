#!/bin/bash

#PBS -N countSig
#PBS -m ea
#PBS -M haoliang.xue@i2bc.paris-saclay.fr
#PBS -q ssfa

dir_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA"

test_set_prefix=$dir_prefix"/b_splitCV"
sig_contig_prefix=$dir_prefix"/d_signatures/d1_merge-first/model.train"
out_dir=$dir_prefix"/d_signatures/d1_merge-first/count-in-tests"

mkdir $out_dir

for n in $(seq 0 4)
do
	/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript /store/USERS/haoliang.xue/development/KaMRaT/related-tools/checkers/contigEval-checker/calc-contig-median.R $sig_contig_prefix$n/signatures-lasso.fa <(zcat $test_set_prefix/kmer-counts.test$n.tsv.gz) $out_dir/contig-counts.test$n.byR.tsv FALSE $test_set_prefix/sampleshuf.test$n.tsv
done

