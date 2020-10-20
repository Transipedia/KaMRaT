#!/bin/bash

usage() { 
	echo -e "Usage: $0 <Required arguments>\nRequired arguments :\n\t-s\tcomplete sample condition file\n\t-p\tpercentage of sample for training\n\t-o\toutput directory\n" 1>&2
	exit
}

[[ $# -eq 0 ]] && usage

while getopts "s:p:o:" opt
do
	case $opt in
		s)	smp_condi_path=$OPTARG
			;;
		p)	train_percent=$OPTARG
			;;
		o)	out_dir=$OPTARG
			;;
		*)	echo "Invalid option provided"
			usage
			;;
	esac
done

if [ "$smp_condi_path" == "" ] || [ "$train_percent" == "" ] || [ "$out_dir" == "" ]
then
	echo
	echo -e "ERROR: required argument is missing"
	echo
	usage
fi

echo
echo "Complete sample condition file:       "$smp_condi_path
echo "Percentage of sample for training:    "$train_percent
echo "Output directory:                     "$out_dir
echo

if ls $out_dir/sampleshuf.tmp.*.tsv 1> /dev/null 2>&1
then
	echo -e "ERROR: some file named sampleshuf.tmp.*.tsv exist in the output directory\n       stop for avoiding conflict"
	exit 1
fi

out_train_path=$out_dir/$(basename $smp_condi_path .tsv).train.tsv
out_test_path=$out_dir/$(basename $smp_condi_path .tsv).test.tsv

# Make sure that sample condtion for TRAIN does not exist
if [ -f $out_dir/sample-condition.train.tsv ]
then
	echo -e "ERROR: sample-condition.train.tsv already exists in the output directory\n       stop for avoiding accidental overwriting"
	exit 1
fi
# Make sure that sample condition for TEST does not exist
if [ -f $out_dir/sample-condition.test.tsv ]
then
	echo -e "ERROR: sample-condition.test.tsv already exists in the output directory\n        stop for avoiding accidental overwriting"
	exit 1
fi

# Splitting starts 
shuf $smp_condi_path | awk -v outd=$out_dir '{print $0 >> outd"/sampleshuf.tmp."$2".tsv"}'

for f in $out_dir/sampleshuf.tmp.*.tsv
do
	row_num=$(wc -l $f | cut -d' ' -f1)
	train_num=`expr $row_num \* $train_percent / 100`
	head -n $train_num $f >> $out_dir/sample-condition.train.tsv
	tail -n `expr $row_num - $train_num` $f >> $out_dir/sample-condition.test.tsv
done

rm -f $out_dir/sampleshuf.tmp.*.tsv
