#!/bin/bash

set -e
pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

tab_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga-relapse/b_splitCV"
out_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA_relapse/c_kamratRes/c1_merge-first"

node_list=(node22 node06 node15 node23 node27)
for i in $(seq 0 4)
do
	out_dir=$out_prefix"/train$i"
	mkdir -p $out_dir

	smp_condi=$tab_dir"/sampleshuf.train"$i".tsv"
	echo "/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat merge -klen 31 -idx-path $out_dir/merge.idx -unstrand -smp-info <(sed '1d' $smp_condi) -interv-method mac:0.32 -quant-mode mean -out-path $out_dir/contig-counts.tsv $tab_dir/kmer-counts.train$i.tsv.gz; rm -f $out_dir/merge.idx" | qsub -N merge-train$i -l nodes=${node_list[$i]}:ppn=1 $pbs_param
done
