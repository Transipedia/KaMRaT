#!/bin/bash

set -e

ensembl_var_1KG="ftp://ftp.ensembl.org/pub/release-101/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz"
out_dir="/data/work/I2BC/haoliang.xue/kamrat/paper/a_ensembl_data"

mkdir $out_dir

cd $out_dir
wget -q $ensembl_var_1KG
gunzip 1000GENOMES-phase_3.gvf.gz

cp /store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa $out_dir/
