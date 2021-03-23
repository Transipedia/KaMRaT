#!/bin/bash

set -e

pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

smp_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga-relapse/b_splitCV"
work_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA_relapse/c_kamratRes/c2_rank-first"

nb_fold=5

node_list=(node22 node13 node15 node23 node14)
for i in $(seq 0 $((nb_fold - 1)))
do
	out_prefix=$work_dir"/train"$i
	smp_condi=$smp_dir"/sampleshuf.train"$i".tsv"
	kmer_count=$smp_dir"/kmer-counts.train"$i".tsv.gz"

	sed "s|\$1|$smp_condi|g; s|\$2|$kmer_count|g; s|\$3|$out_prefix|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga-relapse/d_rank-first/a1_rank_methods_train.qsub | qsub -N rank-train$i -l nodes=${node_list[$i]}:ppn=1 $pbs_param
done
