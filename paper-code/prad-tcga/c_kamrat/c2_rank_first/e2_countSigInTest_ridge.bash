#!/bin/bash

#PBS -N countSig
#PBS -m ea
#PBS -M haoliang.xue@i2bc.paris-saclay.fr
#PBS -q ssfa

dir_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA"

test_set_prefix=$dir_prefix"/b_splitCV"
sig_contig_prefix=$dir_prefix"/d_signatures/d5_rank-first-ridge/model.train"

for n in $(seq 0 4)
do
	/store/USERS/haoliang.xue/development/KaMRaT/bin/contigEvaluate -klen 31 -fasta $sig_contig_prefix$n/signatures-ridge.fa -eval-method median -idx-path $sig_contig_prefix$n/counts.$m.idx -unstrand -smp-info <(sed '1d' $test_set_prefix/sampleshuf.test$n.tsv) -out-path $sig_contig_prefix$n/sig-counts-inTest.median.tsv $test_set_prefix/kmer-counts.test$n.tsv.gz
done

