#!/bin/bash

smp_dir=$1
out_dir=$2
klen=$3

if [ -z "$klen" ]
then
	klen=31
fi

for f1 in $smp_dir/*_1.fasta
do
	smp_name=$(basename $f1 _1.fasta)
	f2=${f1/_1.fasta/_2.fasta}
	/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/jellyfish count -m $klen -s 10000 -t 10 -o $out_dir/$smp_name.jf -F 2 -C $f1 $f2
	/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/jellyfish dump -c $out_dir/$smp_name.jf | sort -k 1 --parallel=10 | /home/haoliang.xue/.conda/envs/dekupl-hlx/bin/pigz -p 10 -c > $out_dir/$smp_name.txt.gz
done

# =======> Jellyfish count
# 	-m, --mer-len=uint32                    *Length of mer 
# 	-s, --size=uint64                       *Initial hash size
# 	-t, --threads=uint32                     Number of threads (1)
# 	-o, --output=string                      Output file (mer_counts.jf)
# 	-F, --Files=uint32                       Number files open simultaneously (1)
#	-C, --canonical                          Count both strand, canonical representation (false)

# =======> Jellyfish dump
# 	-c, --column                             Column format (false)
