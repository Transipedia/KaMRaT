#!/bin/bash

set -e

smp_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga-relapse/b_splitCV"

proj_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA_relapse"
dir_prefix=$proj_dir"/c_kamratRes/c1_merge-first/train"

node_list=(node22 node13 node15 node23 node14)
for i in $(seq 0 4)
do
	smp_condi_path=$smp_dir/sampleshuf.train$i.tsv
	out_dir=$dir_prefix$i
	
	mkdir -p $out_dir

	sed "s|\$1|$smp_condi_path|g; s|\$2|$out_dir|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga-relapse/c_merge-first/c1_selSigs_train.pbs | qsub -N selSig-$i -l nodes=${node_list[$i]}:ppn=4
done

