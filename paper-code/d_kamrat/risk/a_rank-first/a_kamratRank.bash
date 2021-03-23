#!/bin/bash

set -e

pbs_param="-q ssfa -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

smp_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/risk-new"
work_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/ranking-first/risk-new"

nb_fold=5

nlist=(node06 node06)
#wlist=(840565 840566 840567 840568 840569)
for i in $(seq 2 4)
do
	out_prefix=$work_dir"/train"$i
	smp_condi=$smp_dir"/sampleshuf.train"$i".tsv"
	kmer_count=$smp_dir"/kmer-counts.train"$i".tsv.gz"
	for m in `echo snr`
	do
		sed "s|\$1|$smp_condi|g; s|\$2|$kmer_count|g; s|\$3|$m|g; s|\$4|$out_prefix|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/d_kamrat/common/rank_method_train.pbs | qsub -N rsk-rank-$i-$m -l nodes=1:ppn=1 $pbs_param
	done
done
