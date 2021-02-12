#!/bin/bash

read_dir=$1
out_dir=$2
blastdb_name=$3

cat $read_dir/*_1.fasta > $out_dir/all_samples_1.fasta
cat $read_dir/*_2.fasta > $out_dir/all_samples_2.fasta

mkdir -p $out_dir/all_samples

/home/haoliang.xue/.conda/envs/spades/bin/spades.py -o $out_dir/all_samples --rna -1 $out_dir/all_samples_1.fasta -2 $out_dir/all_samples_2.fasta -t 10 -k 31

#	-o <output_dir>             directory to store all the resulting files (required)
#  	--rna                       this flag is required for RNA-Seq data
#	-1 <filename>               file with forward paired-end reads
#	-2 <filename>               file with reverse paired-end reads
#	-t <int>, --threads <int>   number of threads. [default: 16]
#	-k <int> [<int> ...]        list of k-mer sizes (must be odd and less than 128)

echo -e "qseqid\tqlen\tqstart\tqend\tslen\tsstart\tsend\tlength\tnident\tpident\tsseqid" > $out_dir/all_samples/blastn.alignment.tsv
blastn -db $blastdb_name -query $out_dir/all_samples/transcripts.fasta -max_hsps 1 -max_target_seqs 1 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 10 -dust no >> $out_dir/all_samples/blastn.alignment.tsv
