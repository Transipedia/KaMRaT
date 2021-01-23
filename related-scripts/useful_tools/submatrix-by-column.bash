#!/bin/bash

############################################################################################################
# A bash script for retrieving from a given matrix its first column and a given list of other columns      #
# ======================================================================================================== #
# Input:     a matrix to process                                                                           #
#            a file containing a list of column names                                                      #
#            output directory                                                                              #
# Output:    sub-matrix only with the given columns                                                        #
#                                                                                                          #
# Author:    Haoliang Xue                                                                                  #
#            I2BC, CEA, CNRS, Universite Paris-Saclay, Gif sur Yvette, France                              #
############################################################################################################

NOCOLOR='\033[0m'
RED='\033[1;31m'

usage() { 
	echo -e "Usage: $0 <Required arguments>\nRequired arguments :\n\t-m\tthe matrix to process\n\t-l\tfile containing a list of column names\n\t-o\toutput submatrix path\n" 1>&2
	exit
}

[[ $# -eq 0 ]] && usage

while getopts "m:l:o:" opt
do
	case $opt in
		m)	mat_path=$OPTARG
			;;
		l)	col_list_path=$OPTARG
			;;
		o)	out_submat_path=$OPTARG
			;;
		*)	echo -e ${RED}"ERROR: invalid option provided"${NOCOLOR}
			usage
			;;
	esac
done

if [ -z "$mat_path" ] || [ -z "$col_list_path" ] || [ -z "$out_submat_path" ]
then
	echo
	echo -e ${RED}"ERROR: required argument is missing"${NOCOLOR}
	echo
	usage
fi

echo
echo "Origin matrix path:        "$mat_path
echo "Column name list path:     "$col_list_path
echo "Output directory:          "$out_submat_path
echo

col_list_str=`awk '{print $1}' $col_list_path`    # only consider the first column no matter how many columns does the input file have

awk -v col2print="$col_list_str" 'BEGIN {nc2p = split(col2print, c2p);} NR == 1 {for (i = 1; i <= NF; ++i) hname[$i]=i} {printf("%s", $1); for (j = 1; j <= nc2p; ++j) if (c2p[j] in hname) printf("\t%s", $hname[c2p[j]]); print ""}' $mat_path > $out_submat_path
