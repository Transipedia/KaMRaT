#!/bin/bash

set -e

rk_tab_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA/e_gene_level/b_rank_res"
smp_info_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA/b_splitCV/sampleshuf"

node_list=(node22 node06 node15 node23 node27)
for i in $(seq 0 4)
do
	smp_info=$smp_info_prefix.train$i.tsv
	out_dir_prefix=$rk_tab_prefix/train$i-

	sed "s|\$1|$smp_info|g; s|\$2|$out_dir_prefix|g" /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga/e_gene_level/c1_selSignatures.byTrain.qsub | qsub -N selSig-$i -l nodes=${node_list[$i]}:ppn=4
done

