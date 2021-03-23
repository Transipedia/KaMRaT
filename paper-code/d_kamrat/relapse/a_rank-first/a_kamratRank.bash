#!/bin/bash

set -e

pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

smp_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/relapse-new"
work_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/relapse-new"

nb_fold=5
nlist=(node21 node28 node06 node24 node27)
for i in $(seq 0 $((nb_fold - 1)))
do
	out_prefix=$work_dir"/train"$i
	smp_condi=$smp_dir"/sampleshuf.train"$i".tsv"
	kmer_count=$smp_dir"/kmer-counts.train"$i".tsv.gz"
	
	for m in `echo ttest snr lrc nbc` 
	do
		sed "s|\$1|$smp_condi|g; s|\$2|$kmer_count|g; s|\$3|$m|g; s|\$4|$out_prefix|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/d_kamrat/common/rank_method_train.pbs | qsub -N rlps-rank-$i-$m -l nodes=${nlist[$i]}:ppn=1 $pbs_param
	done
done
