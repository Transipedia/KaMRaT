#!/bin/bash

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/g_imoka/a_create_matrix/train"$i

	echo "export IMOKA_MAX_MEM_GB=100; singularity exec -B /store:/store /home/haoliang.xue/tools/iMOKA.img iMOKA_core reduce -i $out_dir/matrix.json -o $out_dir/reduced.matrix" | qsub -N imoka-firstReduce-train$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
