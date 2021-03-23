#!/bin/bash

node_list=(node22 node13 node15 node23 node14)

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA/relapse/train"$i

	echo "export SINGULARITY_BINDPATH=\"/data,/store\"; export IMOKA_MAX_MEM_GB=100; singularity exec /home/haoliang.xue/tools/iMOKA_extended iMOKA_core aggregate -t 70 -T 65 -i $out_dir/reduced.matrix -o $out_dir/aggregated -c $out_dir/matrix.json -m nomap" | qsub -N rlps-agg-train$i -q common -l nodes=${node_list[$i]}:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
