#!/bin/bash

kmer_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga/a_gen_kmer_matrix/kmer_counts"
smp_condi_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga-relapse/b_splitCV/sampleshuf.train"
out_dir_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA_relapse/e_imoka/a_create_matrix/train"

for i in $(seq 0 4)
do
	out_dir=$out_dir_prefix$i
	mkdir -p $out_dir

	awk -v kdir=$kmer_dir 'NR > 1 {print kdir"/"$1".txt\t"$1"\t"$2}' $smp_condi_prefix$i.tsv > $out_dir/create_matrix.tsv
	echo "singularity exec -B /store:/store -B /data:/data /home/haoliang.xue/tools/iMOKA.img iMOKA_core create -i $out_dir/create_matrix.tsv -o $out_dir/matrix.json" | qsub -N rlps-cM-train$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
