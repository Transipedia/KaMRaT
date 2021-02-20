#!/bin/bash

set -e

pbs_param="-q common -m ea -M haoliang.xue@i2bc.paris-saclay.fr"

tab_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA/b_splitCV"
out_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA/c_kamratRes/c2_rank-first"

mkdir -p $out_dir

nb_fold=5

node_list=(node22 node06 node15 node23 node27)
for i in $(seq 0 $((nb_fold - 1)))
do
	smp_condi=$tab_dir"/sampleshuf.train"$i".tsv"
	echo "/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat rank -idx-path $out_dir/rank.train$i.idx -nf-path $out_dir/smp-nf.train$i.tsv -smp-info <(sed '1d' $smp_condi) -score-method snr -out-path $out_dir/ranked-kmer-counts.train$i.tsv $tab_dir/kmer-counts.train$i.tsv.gz; /store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat merge -klen 31 -idx-path $out_dir/merge.train$i.idx -unstrand -smp-info <(sed '1d' $smp_condi) -interv-method mac:0.32 -quant-mode mean -rep-name snr:maxabs -out-path $out_dir/ranked-contig-counts.train$i.tsv <(head -n 100001 $out_dir/ranked-kmer-counts.train$i.tsv)" | qsub -N kRM-train$i -l nodes=${node_list[$i]}:ppn=1 $pbs_param
done
