#!/bin/bash

tab_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/relapse-new"
train_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse-new/train"
out_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse-new/test"

for i in $(seq 3 4)
do
	test_smp_info=$tab_dir/sampleshuf.test$i.tsv
	test_tab=$tab_dir/kmer-counts.test$i.tsv.gz
	
	out_dir=$out_prefix$i
	mkdir -p $out_dir
	
	echo "for m in \`echo ttest snr lrc nbc\`; do in_dir=$train_prefix$i/\$m; /store/USERS/haoliang.xue/development/KaMRaT/bin/contigEvaluate -klen 31 -fasta \$in_dir/signatures-ridge.fa -eval-method median -idx-path $out_dir/counts.\$m.idx -unstrand -smp-info <(sed '1d' $test_smp_info) -out-path $out_dir/contig-counts-ridge.\$m.tsv $test_tab; rm -f $out_dir/counts.\$m.idx; done" | qsub -N rlps-cSig-$i -m ea -M haoliang.xue@i2bc.paris-saclay.fr -q ssfa
done

