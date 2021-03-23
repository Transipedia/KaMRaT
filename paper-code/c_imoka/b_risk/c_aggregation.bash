#!/bin/bash

for i in $(seq 0 4)
do
	out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/iMOKA/risk/train"$i

	echo "export SINGULARITY_BINDPATH=\"/data,/store\"; export IMOKA_MAX_MEM_GB=100; singularity exec /home/haoliang.xue/tools/iMOKA_extended iMOKA_core aggregate -t 70 -T 65 -w 2 -i $out_dir/reduced.matrix -o $out_dir/aggregated -c $out_dir/matrix.json -m nomap" | qsub -N rsk-agg-train$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
