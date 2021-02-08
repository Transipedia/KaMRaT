#!/bin/bash

out_prefix=$1
min_rec=$2
min_rec_abd=$3
smp_list=$4

smp_array=()
for f in "${@:4}"
do
	smp_array+=($f)
done
echo -en "tag\t" > $out_prefix-$min_rec-$min_rec_abd.tsv
echo ${smp_array[@]} | sed 's|/[^ \t]*/sample|sample|g' | sed 's/ /\t/g' | sed 's/.txt.gz//g' >> $out_prefix-$min_rec-$min_rec_abd.tsv
/home/haoliang.xue/.conda/envs/dekupl-hlx/share/dekupl/bin/joinCounts -r $min_rec -a $min_rec_abd ${smp_array[@]} >> $out_prefix-$min_rec-$min_rec_abd.tsv

# =======> DE-kupl_joinCounts
#	-r INT    min recurrence [1]
# 	-a INT    min recurrence abundance [1]	
