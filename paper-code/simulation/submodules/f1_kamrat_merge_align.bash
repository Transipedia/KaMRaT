#!/bin/bash

in_tab_path=$1
out_prefix=$2
interv_method=$3
blastdb_name=$4
klen=$5
min_overlap=$6

if [ ! -z "$min_overlap" ]
then
	min_overlap="-min-overlap $min_overlap"
fi

in_tab_name=$(basename $in_tab_path .tsv)
	
out_tab_path=$out_prefix-$in_tab_name.${interv_method/:/_}.tsv
out_idx_path=$out_prefix-$in_tab_name.${interv_method/:/_}.count.idx
/store/USERS/haoliang.xue/development/KaMRaT/bin/kamrat merge -klen $klen -idx-path $out_idx_path $min_overlap -unstrand -interv-method $interv_method -quant-mode mean -out-path $out_tab_path $in_tab_path

out_fa_path=$out_prefix-${in_tab_name/raw-counts/contigs}.${interv_method/:/_}.fa
awk 'NR > 1 {print ">contig_"NR - 1"\n"$1}' $out_tab_path > $out_fa_path

out_align_path=${out_tab_path/raw-counts/contigs-align}
echo -e "qseqid\tqlen\tqstart\tqend\tslen\tsstart\tsend\tlength\tnident\tpident\tsseqid" > $out_align_path
blastn -db $blastdb_name -query $out_fa_path -max_hsps 1 -max_target_seqs 1 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 10 -dust no >> $out_align_path
