#!/bin/bash

node_list=(node22 node13 node15 node23 node14)

work_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/g_imoka"

for i in $(seq 0 4)
do
	mat_dir=$work_dir"/a_create_matrix/train"$i
	out_dir=$work_dir"/b_rf_models/train"$i

	mkdir -p $out_dir

	echo "singularity exec -B /store:/store -B /data:/data /home/haoliang.xue/tools/iMOKA.img random_forest.py -m 100 -t 8 $mat_dir/aggregated.kmers.matrix $out_dir/randomforest" | qsub -N imoka-rf-train$i -q common -l nodes=${node_list[$i]}:ppn=8 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
