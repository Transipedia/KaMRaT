#!/bin/bash

############################################################################################################
# A bash script for splitting a sample-condition file into $n$ training-testing pairs for cross-validation #
# ======================================================================================================== #
# Input:     a sample-condition file without header line                                                   #
#            number of fold for cross-validation ($n$)                                                     #
#            output directory                                                                              #
# Output:    $n$ training-testing pairs of sample-condition files                                          #
#                sampleshuf.train$n$.tsv is associated with sampleshuf.test$n$.tsv                         #
#                                                                                                          #
# For avoiding removing files by accident,                                                                 #
# please make sure that the output directory is empty before running the script                            #
#                                                                                                          #
# Author:    Haoliang Xue                                                                                  #
#            I2BC, CEA, CNRS, Universite Paris-Saclay, Gif sur Yvette, France                              #
############################################################################################################

usage() { 
	echo -e "Usage: $0 <Required arguments>\nRequired arguments :\n\t-s\tcomplete sample condition file\n\t-n\tfold number for cross-validation\n\t-o\toutput directory\n" 1>&2
	exit
}

[[ $# -eq 0 ]] && usage

while getopts "s:n:o:" opt
do
	case $opt in
		s)	smp_condi_path=$OPTARG
			;;
		n)	nb_fold=$OPTARG
			;;
		o)	out_dir=$OPTARG
			;;
		*)	echo "Invalid option provided"
			usage
			;;
	esac
done

if [ -z "$smp_condi_path" ] || [ -z "$nb_fold" ] || [ -z "$out_dir" ]
then
	echo
	echo -e "ERROR: required argument is missing"
	echo
	usage
fi

if [[ ! $nb_fold =~ ^[1-9][0-9]*$ ]] || [ $nb_fold -eq 1 ]
then
	echo
	echo -e "ERROR: invalid fold number: "$nb_fold
	echo
	exit
fi

echo
echo "Complete sample condition file:       "$smp_condi_path
echo "Fold number for cross-validation:     "$nb_fold
echo "Output directory:                     "$out_dir
echo

get_seeded_random() {
	seed="$1"
	openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

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
shuf --random-source=<(get_seeded_random 91400) $smp_condi_path | awk -v outd=$out_dir -v nfold=$nb_fold '{print $0 >> outd"/sampleshuf.tmp."$2"_"NR%nfold".tsv"}'

for i in $(seq 0 $(( $nb_fold - 1 )))
do
	cat $out_dir/sampleshuf.tmp.*_$i.tsv > $out_dir/sampleshuf.tmp.fold$i.tsv
done

for i in $(seq 0 $(( $nb_fold - 1 )))
do
	cat `ls $out_dir/sampleshuf.tmp.fold*.tsv | grep -v sampleshuf.tmp.fold$i.tsv` > $out_dir/sampleshuf.train$i.tsv
	cat $out_dir/sampleshuf.tmp.fold$i.tsv > $out_dir/sampleshuf.test$i.tsv
done

rm -f $out_dir/sampleshuf.tmp.*.tsv
