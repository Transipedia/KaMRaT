#!/bin/bash

set -e
pbs_param="-q ssfa -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

tab_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/relapse-new"
out_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse-new"

for i in $(seq 0 4)
do
	out_dir=$out_prefix"/train$i"
	mkdir -p $out_dir

	smp_condi=$tab_dir"/sampleshuf.train"$i".tsv"
	echo "/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat merge -klen 31 -idx-path $out_dir/merge.idx -unstrand -smp-info <(sed '1d' $smp_condi) -interv-method mac:0.32 -quant-mode mean -out-path $out_dir/merged-counts.tsv $tab_dir/kmer-counts.train$i.tsv.gz; rm -f $out_dir/merge.idx" | qsub -N rlps-merge-$i -l nodes=1:ppn=1 $pbs_param
done
