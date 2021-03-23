#!/bin/bash

set -e
pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr -W depend=afterok:840718"

tab_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/risk-new"
out_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/risk-new"

for i in $(seq 3 3)
do
	out_dir=$out_prefix"/train$i"
	mkdir -p $out_dir

	smp_condi=$tab_dir"/sampleshuf.train"$i".tsv"
	echo "/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat merge -klen 31 -idx-path $out_dir/merge.idx -unstrand -smp-info <(sed '1d' $smp_condi) -interv-method mac:0.32 -quant-mode mean -out-path $out_dir/merged-counts.tsv $tab_dir/kmer-counts.train$i.tsv.gz; rm -f $out_dir/merge.idx" | qsub -N rsk-merge-$i -l nodes=node24:ppn=1 $pbs_param
done
