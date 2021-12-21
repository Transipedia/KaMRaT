#!/bin/bash

countTab=$1

awk 'BEGIN {count0 = 0} NR > 1 {for (i = 2; i <= NF; i ++) if ($i == 0) count0 ++;} END {print "Total count number: "(NF - 1) * (NR - 1); print "Zero number: "count0;}' $countTab

