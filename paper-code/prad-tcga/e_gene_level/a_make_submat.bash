#!/bin/bash

#PBS -N make_submat
#PBS -m ea
#PBS -M haoliang.xue@i2bc.paris-saclay.fr

init_gene_tab="/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Prostate_TCGA_LRvsHR/DEkupl_LRvsHR_limma_Mask/Result_gene_counts/gene_counts.tsv.gz"
smp_info_prefix="/data/work/I2BC/haoliang.xue/PRAD_TCGA/b_splitCV/sampleshuf"
out_dir="/data/work/I2BC/haoliang.xue/PRAD_TCGA/e_gene_level"

mkdir $out_dir

for i in $(seq 0 4)
do
	/store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/selectSubmat.bash -m <(zcat $init_gene_tab) -l <(sed '1d' $smp_info_prefix.train$i.tsv) -o $out_dir/gene-counts.train$i.tsv
	/store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/selectSubmat.bash -m <(zcat $init_gene_tab) -l <(sed '1d' $smp_info_prefix.test$i.tsv) -o $out_dir/gene-counts.test$i.tsv
done
