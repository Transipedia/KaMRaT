#!/bin/bash

proj_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/kamrat-paper/prad-tcga"
smp_prefix=$proj_dir"/h_mergef-deseq2/data/smp_condi_full.train"
out_dir_prefix=$proj_dir"/h_mergef-deseq2/train"

pbs_param="-m ea -M haoliang.xue@i2bc.paris-saclay.fr -q common"

nodelist=(node14 node15 node22 node23 node27)

echo ${nodelist[@]}

for i in $(seq 0 4)
do
	s=$smp_prefix$i".tsv"
	out_dir=$out_dir_prefix$i
	mkdir -p $out_dir"/tmp"
	f=$out_dir_prefix$i"/contig-counts.tsv.gz"
	echo "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript /home/haoliang.xue/.conda/envs/dekupl-hlx/share/dekupl/bin/DESeq2_diff_method.R /home/haoliang.xue/.conda/envs/dekupl-hlx/share/dekupl/bin/TtestFilter $f $s 1 0 LR HR 8 1000000 $out_dir/tmp/ $out_dir/diff-counts.tsv.gz $out_dir/raw_pvals.txt.gz $out_dir/tmp/log fixed" | qsub -N deseq2-train$i $pbs_param -l nodes=${nodelist[$i]}:ppn=8
	echo ""
done

