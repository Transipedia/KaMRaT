#!/bin/bash

kamrat_tab_path=$1
out_tab_path=$2

awk 'NR == 1 {printf("%s\t%s\t%s\tscore\tmeanA\tmeanB\tlog2FC", $2, $1, $3); for (i = 9; i <= NF; i++) printf("\t%s", $i); print ""} NR > 1 {printf("%s\t%s\t%s\t%s\t%s\t%s", $2, $1, $3, $4, $5, $6); printf("\t%.6f", log(($6 + 1)/($5 + 1))/log(2)); for (i = 9; i <= NF; i++) printf("\t%s", $i); print ""}' $kamrat_tab_path | sed 's/[Ii]nf/999/g' > $out_tab_path
