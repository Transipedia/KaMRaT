#!/bin/bash

node_list=(node22 node13 node15 node23 node14)

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/g_imoka/a_create_matrix/train"$i

	echo "export IMOKA_MAX_MEM_GB=100; singularity exec -B /store:/store -B /data:/data /home/haoliang.xue/tools/iMOKA.img iMOKA_core aggregate -i $out_dir/reduced.matrix -o $out_dir/aggregated -c $out_dir/matrix.json -m default" | qsub -N imoka-agg-train$i -q common -l nodes=${node_list[$i]}:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
