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

NOCOLOR='\033[0m'
RED='\033[1;31m'

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
		*)	echo -e ${RED}"ERROR: invalid option provided"${NOCOLOR}
			usage
			;;
	esac
done

if [ -z "$smp_condi_path" ] || [ -z "$nb_fold" ] || [ -z "$out_dir" ]
then
	echo
	echo -e ${RED}"ERROR: required argument is missing"${NOCOLOR}
	echo
	usage
fi

if [[ ! $nb_fold =~ ^[1-9][0-9]*$ ]] || [ $nb_fold -eq 1 ]
then
	echo
	echo -e ${RED}"ERROR: invalid fold number: "$nb_fold${NOCOLOR}
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

if ls $out_dir/sampleshuf.tmp*.tsv 1> /dev/null 2>&1
then
	echo -e ${RED}"ERROR: some file named sampleshuf.tmp*.tsv exist in the output directory\n       stop for avoiding conflict"${NOCOLOR}
	exit 1
fi

out_train_path=$out_dir/$(basename $smp_condi_path .tsv).train.tsv
out_test_path=$out_dir/$(basename $smp_condi_path .tsv).test.tsv

# Make sure that sample condtion for TRAIN does not exist
if [ -f $out_dir/sampleshuf.train*.tsv ]
then
	echo -e ${RED}"ERROR: sampleshuf.train*.tsv already exists in the output directory\n       stop for avoiding accidental overwriting"${NOCOLOR}
	exit 1
fi
# Make sure that sample condition for TEST does not exist
if [ -f $out_dir/sampleshuf.test*.tsv ]
then
	echo -e ${RED}"ERROR: sampleshuf.test*.tsv already exists in the output directory\n        stop for avoiding accidental overwriting"${NOCOLOR}
	exit 1
fi

# Splitting starts
awk -v outd=$out_dir '{print $0 >> outd"/sampleshuf.tmp1."$2".tsv"}' $smp_condi_path
for fs in $out_dir/sampleshuf.tmp1.*.tsv
do
	shuf --random-source=<(get_seeded_random 91400) $fs | awk -v outd=$out_dir -v nfold=$nb_fold '{print $0 >> outd"/sampleshuf.tmp2."$2"_"NR%nfold".tsv"}'
done

for i in $(seq 0 $(( $nb_fold - 1 )))
do
	cat $out_dir/sampleshuf.tmp2.*_$i.tsv > $out_dir/sampleshuf.tmp3.fold$i.tsv
done

for i in $(seq 0 $(( $nb_fold - 1 )))
do
	cat `ls $out_dir/sampleshuf.tmp3.fold*.tsv | grep -v sampleshuf.tmp3.fold$i.tsv` > $out_dir/sampleshuf.train$i.tsv
	cat $out_dir/sampleshuf.tmp3.fold$i.tsv > $out_dir/sampleshuf.test$i.tsv
done

rm -f $out_dir/sampleshuf.tmp*.tsv
