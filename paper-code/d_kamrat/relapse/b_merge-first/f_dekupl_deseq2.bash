#!/bin/bash

dkpl_smp_condi="/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Prostate_TCGA_relapse/DEkupl_run_TCGA_relapse/metadata/sample_conditions_full.tsv"
smp_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/CVMatrices/relapse-new/sampleshuf.train"
out_dir_prefix="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/PRAD_TCGA/merging-first/relapse-new/train"

pbs_param="-m ea -M haoliang.xue@i2bc.paris-saclay.fr -q common"

for i in $(seq 0 4)
do
	f=$out_dir_prefix$i"/merged-counts.tsv"
	out_dir=$out_dir_prefix$i"/deseq2"
	mkdir -p $out_dir"/tmp"
	
	slist=`awk 'NR > 1 {print $1}' $smp_prefix$i.tsv`
	awk -v smplist="$slist" 'BEGIN {ns = split(smplist, sl)} NR == 1 {print $0} NR > 1 {rowdict[$1]=$0} END {for (i = 1; i <= ns; i++) print rowdict[sl[i]]}' $dkpl_smp_condi > $out_dir/sample_condition_full.tsv 
	
	echo "awk 'NR == 1 {printf(\"%s\", \$3); for (i = 4; i <= NF; i++) printf(\"\t%s\", \$i); print \"\"} NR > 1 {printf(\"%s\", \$3); for (i = 4; i <= NF; i++) printf(\"\t%d\", \$i + 0.5); print \"\"}' $f | /home/haoliang.xue/.conda/envs/dekupl-hlx/bin/pigz -p 8 -c > $out_dir/merged-counts.tsv.gz; /home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript /home/haoliang.xue/.conda/envs/dekupl-hlx/share/dekupl/bin/DESeq2_diff_method.R /home/haoliang.xue/.conda/envs/dekupl-hlx/share/dekupl/bin/TtestFilter $out_dir/merged-counts.tsv.gz $out_dir/sample_condition_full.tsv 1 0 No Yes 8 1000000 $out_dir/tmp/ $out_dir/diff-counts.tsv.gz $out_dir/raw_pvals.txt.gz $out_dir/tmp/log fixed" | qsub -N rlps-deseq2-train$i $pbs_param -l nodes=1:ppn=8
done

