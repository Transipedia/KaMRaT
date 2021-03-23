#!/bin/bash

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/g_imoka/a_create_matrix/train"$i
	mkdir -p $out_dir

	echo "singularity exec -B /store:/store /home/haoliang.xue/tools/iMOKA.img iMOKA_core create -i /store/USERS/haoliang.xue/development/KaMRaT/paper-code/prad-tcga/g_imoka/data/create_matrix.train$i.tsv -o $out_dir/matrix.json" | qsub -N imoka-createMat-train$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
