#!/bin/bash

set -e

dir_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA_relapse"

tab_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga-relapse/b_splitCV"
out_prefix=$dir_prefix"/c_kamratRes/c1_merge-first"

node_list=(node22 node13 node15 node23 node27)
for i in $(seq 0 4)
do
	sed "s|\$1|$tab_dir/kmer-counts.test$i.tsv.gz|g; s|\$2|$tab_dir/sampleshuf.test$i.tsv|g; s|\$3|$out_prefix/train$i|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga-relapse/c_merge-first/d1_count_test_train.qsub | qsub -N countTest-$i -l nodes=${node_list[$i]}:ppn=1
done

