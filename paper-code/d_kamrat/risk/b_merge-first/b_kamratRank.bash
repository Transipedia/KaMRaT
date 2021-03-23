#!/bin/bash

set -e

pbs_param="-q ssfa -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

smp_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/risk-new"
work_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/risk-new/train"

nb_fold=5

jlist=(840597 840598 840599 840600 840601)
for i in $(seq 0 $((nb_fold - 1)))
do
	smp_condi=$smp_dir"/sampleshuf.train"$i".tsv"
	in_tab=$work_prefix$i"/merged-counts.tsv"
	out_prefix=$work_prefix$i
	for m in `echo ttest snr lrc nbc`
	do
		sed "s|\$1|$smp_condi|g; s|\$2|$in_tab|g; s|\$3|$m|g; s|\$4|$out_prefix|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/d_kamrat/common/rank_method_train.pbs | qsub -N rsk-rank-$i-$m -l nodes=1:ppn=1 -W depend=afterok:${jlist[$i]} $pbs_param
	done
done
