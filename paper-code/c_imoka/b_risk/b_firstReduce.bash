#!/bin/bash

node_list=(node22 node13 node15 node23 node14)

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA/risk/train"$i

	echo "export OMP_NUM_THREADS=8; export IMOKA_MAX_MEM_GB=100; singularity exec -B /store:/store -B /data:/data /home/haoliang.xue/tools/iMOKA_extended iMOKA_core reduce -i $out_dir/matrix.json -o $out_dir/reduced.matrix" | qsub -N rsk-fR-train$i -q common -l nodes=${node_list[$i]}:ppn=8 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
