#/bin/bash

#PBS -N gen_events
#PBS -q common
#PBS -l nodes=1:ppn=1
#PBS -m bea
#PBS -M haoliang.xue@i2bc.paris-saclay.fr

gvf_path=$1
# gvf_path=/home/haoliang.xue/media/data/kamrat/paper/a_ensembl_data/1000GENOMES-phase_3.gvf
out_dir=$2
#out_dir=/home/haoliang.xue/media/data/kamrat/paper/b_reference_gen
rand_seed=91400

if ls $out_dir/b_tmp.1000GENOMES-phase_3.*.gvf 1> /dev/null 2>&1
then
	echo "ERROR: tmp file exists"
	exit 1
fi

get_seeded_random() {
	seed="$1"
	openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

# Divide GVF file into categories by event type
grep ^[^#] $gvf_path | awk -v awk_out_dir=$out_dir '{print $0 >> awk_out_dir"/b_tmp.1000GENOMES-phase_3."$3".gvf"}'

# Shuffle 100 variations for each event category
for f in $out_dir/b_tmp.1000GENOMES-phase_3.{SNV,deletion,insertion,indel}.gvf; do
	shuf -n100 --random-source=<(get_seeded_random $rand_seed) $f >> $out_dir/b_1000GENOMES-phase_3.sel.gvf
done

rm -f $out_dir/b_tmp.1000GENOMES-phase_3.*.gvf
