#!/bin/bash

node_list=(node22 node13 node15 node23 node14)

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA/relapse/train"$i

	echo "export SINGULARITY_BINDPATH=\"/data,/store\"; singularity exec /home/haoliang.xue/tools/iMOKA_extended random_forest.py -m 5000 $out_dir/aggregated.kmers.matrix $out_dir/randomforest" | qsub -N rlps-rf-train$i -q common -l nodes=${node_list[$i]}:ppn=2 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
