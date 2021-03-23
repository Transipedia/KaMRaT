#!/bin/bash


full_smp_condi_path="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/data/sample-condition.risk.tsv"
jellyfish_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/data/jellyfish_res/kmer_counts"
out_dir="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/risk-new"

nb_fold=5

mkdir -p $out_dir

smp_condi_header=`head -n1 $full_smp_condi_path`
/store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/splitCV.bash -s <(sed '1d' $full_smp_condi_path) -n $nb_fold -o $out_dir

for i in $(seq 0 $((nb_fold - 1)))
do
	f=$out_dir/sampleshuf.train$i.tsv
	echo "set -e; t=${f/sampleshuf/kmer-counts}.gz; awk 'BEGIN{printf(\"tag\")} {printf(\"\t%s\", \$1)} END{print \"\"}' $f | gzip -c >> \$t; path_list=\`awk -v jd=$jellyfish_dir '{print jd\"/\"\$1\".txt.gz\"}' $f\`; /store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/dekupl-joinCounts/joinCounts -r 1 -a 5 \${path_list[@]} | gzip -c >> \$t; sed -i \"1i $smp_condi_header\" $f" | qsub -N makeCV-risk-train$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done

for i in $(seq 0 $((nb_fold - 1)))
do
	f=$out_dir/sampleshuf.test$i.tsv
	echo "set -e; t=${f/sampleshuf/kmer-counts}.gz; awk 'BEGIN{printf(\"tag\")} {printf(\"\t%s\", \$1)} END{print \"\"}' $f | gzip -c >> \$t; path_list=\`awk -v jd=$jellyfish_dir '{print jd\"/\"\$1\".txt.gz\"}' $f\`; /store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/dekupl-joinCounts/joinCounts -r 1 -a 1 \${path_list[@]} | gzip -c >> \$t; sed -i \"1i $smp_condi_header\" $f" | qsub -N makeCV-risk-test$i -q ssfa -l nodes=1:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr
done
