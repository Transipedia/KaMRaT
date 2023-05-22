#!/bin/bash
# Author: Haoliang Xue

set -e

# Involved programs
KAMRAT_IMG="/home/haoliang.xue_ext/KaMRaT.sif" # need to change the singularity image path accordingly
TIME='/usr/bin/time -f "Command - %C\n\tuser time in seconds -  %U\n\tsystem (kernel) time in seconds - %S\n\telapsed real time (wall clock) in seconds - %e\n\tpercent of CPU this job got - %P"'

# Inputs
IN_DIR="/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/LUADseo/"
IN_MAT=$IN_DIR"kmer-counts.allsmp.tsv.gz"
IN_META=$IN_DIR"all-samples.tsv"
IN_NFFILE=$IN_DIR"nf-values.allsmp.tsv"

# Input examples
# IN_DIR="/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/toyroom/data/"
# IN_MAT=$IN_DIR"kmer-counts.subset4toy.tsv.gz"
# IN_META=$IN_DIR"sample-condition.toy.tsv"
# IN_NFFILE=$IN_DIR"nf_values.txt"

# Outputs
OUT_DIR="/data/work/I2BC/haoliang.xue/kamrat-new-res/test-bashscript/" # need to change accordingly
IDX_DIR=$OUT_DIR"index/"

mkdir -p $IDX_DIR

bash -c "$TIME apptainer exec -B '/store:/store' -B '/data:/data' $KAMRAT_IMG kamrat index -intab $IN_MAT -outdir $IDX_DIR -klen 31 -unstrand -nffile $IN_NFFILE &> $IDX_DIR/log-index.txt"

bash -c "$TIME apptainer exec -B '/store:/store' -B '/data:/data' $KAMRAT_IMG kamrat merge -idxdir $IDX_DIR -overlap 30-15 -interv pearson:0.20 -outpath $OUT_DIR/merged-res.bin &> $OUT_DIR/log-merge-rank.txt" 

bash -c "$TIME apptainer exec -B '/store:/store' -B '/data:/data' $KAMRAT_IMG kamrat rank -idxdir $IDX_DIR -rankby ttest.padj -design $IN_META -outpath $OUT_DIR/merge-rank-res.tsv -with $OUT_DIR/merged-res.bin:mean -withcounts &>> $OUT_DIR/log-merge-rank.txt"

