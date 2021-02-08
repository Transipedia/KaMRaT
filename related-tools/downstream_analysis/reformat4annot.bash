#!/bin/bash

#########################################################################
# A bash script to reformat KaMRaT-merge results for DE-kupl annotation #
# --------------------------------------------------------------------- #
# note: only binary condition acceptable (requires log2fc) !            #
#                                                                       #
# author: Haoliang Xue (I2BC)                                           #
#########################################################################

kamrat_tab_path=$1
smp_info_path=$2
out_tab_path=$3

awk 'NR == 1 {printf("%s\t%s\t%s\tscore\tmeanA\tmeanB\tlog2FC", $2, $1, $3); for (i = 9; i <= NF; i++) printf("\t%s", $i); print ""} NR > 1 {printf("%s\t%s\t%s\t%s\t%s\t%s", $2, $1, $3, $4, $5, $6); printf("\t%.6f", log(($6 + 1)/($5 + 1))/log(2)); for (i = 9; i <= NF; i++) printf("\t%s", $i); print ""}' $kamrat_tab_path > $out_tab_path
