#!/bin/bash

kamrat_tab_path=$1
out_tab_path=$2

awk 'NR == 1 {printf("%s\tpvalue\t%s\t%s\tlog2FC", $1, $(NF - 1), $NF); for (i = 3; i <= NF - 2; i++) printf("\t%s", $i); print ""} NR > 1 {printf("%s\t%f\t%f\t%f\t%f", $1, $2, $(NF - 1), $NF, log($NF/$(NF - 1))/log(2)); for (i = 3; i <= NF - 2; i++) printf("\t%s", $i); print ""}' $kamrat_tab_path > $out_tab_path
