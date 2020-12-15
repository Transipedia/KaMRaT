#!/bin/bash

usage() {
	echo -e "Usage:\n\tbash $(basename $0) columns-to-select count-table-path" 1>&2
	exit
}

[[ $# -eq 0 ]] && usage

sub_smp_arr=($(awk '{print $1}' $1))
count_tab_path=$2

sub_smp_str=" "${sub_smp_arr[*]}" "

awk -v subsmp="$sub_smp_str" '{if (NR == 1) split($0, colnames); printf("%s\t", $1); for (i = 1; i <= length(colnames); ++i) if (match(subsmp, colnames[i])) printf("\t%s", $i); print ""}' $count_tab_path
