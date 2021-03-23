#!/bin/bash

set -e

proj_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga"

nodelist=(node14 node15 node22 node23 node27)
echo ${nodelist[@]}

for i in $(seq 0 4)
do
	smp_condi=$proj_dir"/b_splitCV/sampleshuf.train"$i".tsv"
	merge_tab=$proj_dir"/c_kamratRes/c1_merge-first/contig-counts.train$i.tsv"
	out_dir=$proj_dir"/h_mergef-deseq2/train"$i""
	echo "/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat rank -idx-path $out_dir/rank.snr.idx -nf-path $out_dir/smp-nf.snr.tsv -smp-info <(sed '1d' $smp_condi) -score-method snr -out-path $out_dir/snr-ranked-contig-counts.tsv $merge_tab" | qsub -N kamratSNR.train$i -m ea -M haoliang.xue@i2bc.paris-saclay.fr -q common -l nodes=${nodelist[$i]}:ppn=1
done

