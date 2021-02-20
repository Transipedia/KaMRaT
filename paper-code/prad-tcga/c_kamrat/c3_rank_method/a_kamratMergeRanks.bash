#!/bin/bash

set -e

pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

smp_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA/b_splitCV"
tab_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA/c_kamratRes/c1_merge-first"
out_dir_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA/c_kamratRes/c3_rank-methods"

nb_fold=5

node_list=(node22 node06 node15 node23 node27)
for i in $(seq 0 $((nb_fold - 1)))
do
	out_dir=$out_dir_prefix"/train"$i
	mkdir -p $out_dir
	smp_condi=$smp_dir"/sampleshuf.train"$i".tsv"
	ctg_count=$tab_dir"/contig-counts.train"$i".tsv"

	sed "s|\$1|$smp_condi|g; s|\$2|$ctg_count|g; s|\$3|$out_dir|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga/c_kamrat/c3_rank_method/a1_rank_per_train.pbs | qsub -N kR-train$i -l nodes=${node_list[$i]}:ppn=1 $pbs_param
done
