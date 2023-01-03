#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile c_rnaspades.smk --cluster "qsub -q common -l nodes=node22:ppn=6 -l mem=120G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 1 -p --latency-wait 60 --rerun-incomplete >> workflow_c_rnaspades.txt &

# Involved programs
SPADES = "/home/haoliang.xue_ext/SPAdes-3.15.5-Linux/bin/spades.py"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"

# Inputs
BLAST_DB = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Data/bench-merge/gencode.v34.transcripts.fa"

# Outputs
RES_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
PLSTR_DIR = RES_DIR + "polyester_res/"
SPADES_DIR = RES_DIR + "spades_res/"

SMP_LIST = ["sample_01", "sample_02", "sample_03", "sample_04", "sample_05", "sample_06", "sample_07", "sample_08", "sample_09", "sample_10", "sample_11", "sample_12", "sample_13", "sample_14", "sample_15", "sample_16", "sample_17", "sample_18", "sample_19", "sample_20"]
MEAN_DEPTH_LIST = [0.1, 0.5, 1, 2, 5, 10]

# ===== Workflow ===== #
rule all:
    input:
        expand(SPADES_DIR + "depth_{meandepth}/allreads/blastn_align.tsv", meandepth = MEAN_DEPTH_LIST),
        expand(SPADES_DIR + "depth_{meandepth}/allkmers/blastn_align.tsv", meandepth = MEAN_DEPTH_LIST)

rule spades_allreads:
    input:
        fa1 = expand(PLSTR_DIR + "depth_{meandepth}/{smp}_1.fasta", meandepth = ["{meandepth}"], smp = SMP_LIST),
        fa2 = expand(PLSTR_DIR + "depth_{meandepth}/{smp}_2.fasta", meandepth = ["{meandepth}"], smp = SMP_LIST)
    output:
        SPADES_DIR + "depth_{meandepth}/allreads/transcripts.fasta"
    log:
        SPADES_DIR + "depth_{meandepth}/allreads/log-spades.txt"
    threads: 6
    params:
        fa1 = SPADES_DIR + "depth_{meandepth}/allreads_1.fasta",
        fa2 = SPADES_DIR + "depth_{meandepth}/allreads_2.fasta",
        spades_out = SPADES_DIR + "depth_{meandepth}/allreads/"
    shell:
        """
        cat {input.fa1} > {params.fa1}
        cat {input.fa2} > {params.fa2}
        {SPADES} -o {params.spades_out} --rna -1 {params.fa1} -2 {params.fa2} -t 6 &> {log}
        """

rule spades_allkmers:
    input:
        RES_DIR + "matrices/depth_{meandepth}/kmer-counts.tsv"
    output:
        SPADES_DIR + "depth_{meandepth}/allkmers/transcripts.fasta"
    log:
        SPADES_DIR + "depth_{meandepth}/allkmers/log-spades.txt"
    threads: 6
    params:
        fa = SPADES_DIR + "depth_{meandepth}/allkmers.fasta",
        spades_out = SPADES_DIR + "depth_{meandepth}/allkmers/"
    shell:
        """
        awk 'NR > 1 {{print ">kmer_"NR - 1"\\n"$1}}' {input} > {params.fa}
        {SPADES} -o {params.spades_out} --rna --s 1 {params.fa} -t 6 &> {log}
        """

rule blastn:
    input:
        SPADES_DIR + "depth_{meandepth}/{mode}/transcripts.fasta"
    output:
        SPADES_DIR + "depth_{meandepth}/{mode}/blastn_align.tsv"
    threads: 6
    log:
        SPADES_DIR + "depth_{meandepth}/{mode}/log-blastn.txt"
    shell:
        """
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {BLAST_DB} -query {input} -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 6 -dust no >> {output} 2> {log}
        """


# ============ SPAdes parameters ============ #
#	-o <output_dir>             directory to store all the resulting files (required)
#  	--rna                       this flag is required for RNA-Seq data
#	-1 <filename>               file with forward paired-end reads
#	-2 <filename>               file with reverse paired-end reads
#       --s <#> <filename>          file with unpaired reads for single reads library number <#>.
#	-t <int>, --threads <int>   number of threads. [default: 16]
#	-k <int> [<int> ...]        list of k-mer sizes (must be odd and less than 128)
