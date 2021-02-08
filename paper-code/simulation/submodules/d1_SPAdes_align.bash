#!/bin/bash

read_dir=$1
out_dir=$2
blastdb_name=$3

for f1 in $read_dir/*_1.fasta
do
	smp=$(basename $f1 _1.fasta)
	f2=${f1/_1.fasta/_2.fasta}
	mkdir -p $out_dir/$smp
	/home/haoliang.xue/.conda/envs/spades/bin/spades.py -o $out_dir/$smp/ --rna -1 $f1 -2 $f2 -t 10 -k 31
done

#	-o <output_dir>             directory to store all the resulting files (required)
#  	--rna                       this flag is required for RNA-Seq data
#	-1 <filename>               file with forward paired-end reads
#	-2 <filename>               file with reverse paired-end reads
#	-t <int>, --threads <int>   number of threads. [default: 16]
#	-k <int> [<int> ...]        list of k-mer sizes (must be odd and less than 128)

for f in $out_dir/sample_[0-9][0-9]/transcripts.fasta
do
	smp_name=$(basename ${f%/*})
	smp_dir=$(dirname $f)
	echo -e "qseqid\tqlen\tqstart\tqend\tslen\tsstart\tsend\tlength\tnident\tpident\tsseqid" > $smp_dir/blastn.alignment.tsv
	blastn -db $blastdb_name -query $f -max_hsps 1 -max_target_seqs 1 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 10 -dust no >> $smp_dir/blastn.alignment.tsv
done
