#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile check_b_kmer_error.smk --cluster "qsub -q common -l nodes=1:ppn=6 -l mem=20G -l walltime=300:00:00" --jobs 8 -p --latency-wait 60 --rerun-incomplete >> workflow_check_b_kmer_error.txt &


# Parameters
MEAN_DEPTH_LIST = [1, 5, 10, 20, 30, 40, 50]

# Involved programs
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"

# Inputs and outputs
REF_FASTA = "/store/plateformes/CALCUL/SSFA_KaMRaT/Data/gc34.50-53.fa"
RES_DIR = f"/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
BLAST_DB = RES_DIR + "blast_db/"
MAT_DIR = RES_DIR + "matrices/"
CHECK_KMER_DIR = RES_DIR + "check_kmer_error/"


# ===== Workflow ===== #
rule all:
    input:
        expand(CHECK_KMER_DIR + "err-free_1-1/depth_{meandepth}.align.tsv",
               meandepth = MEAN_DEPTH_LIST),
        expand(CHECK_KMER_DIR + "err-illumina5_2-1/depth_{meandepth}.align.tsv",
               meandepth = MEAN_DEPTH_LIST),
        expand(CHECK_KMER_DIR + "err-illumina5_2-5/depth_{meandepth}.align.tsv",
               meandepth = MEAN_DEPTH_LIST)

rule makeblastdb:
    input:
        REF_FASTA
    output:
        BLAST_DB + "gc34.50-53.fa",
        BLAST_DB + "gc34.50-53.fa.nhr",
        BLAST_DB + "gc34.50-53.fa.nin",
        BLAST_DB + "gc34.50-53.fa.nog",
        BLAST_DB + "gc34.50-53.fa.nsd",
        BLAST_DB + "gc34.50-53.fa.nsi",
        BLAST_DB + "gc34.50-53.fa.nsq"
    threads: 1
    shell:
        """
        cp {REF_FASTA} {BLAST_DB}
        {MKBLASTDB} -in {input} -dbtype nucl -parse_seqids -out {BLAST_DB}gc34.50-53.fa
        """

rule blastn:
    input:
        kmer_mat = MAT_DIR + "{mode}/depth_{meandepth}/kmer-counts-{min_rec}-{min_abd}.tsv",
        blast_db = BLAST_DB + "gc34.50-53.fa"
    output:
        CHECK_KMER_DIR + "{mode}_{min_rec}-{min_abd}/depth_{meandepth}.align.tsv"
    params:
        CHECK_KMER_DIR + "{mode}_{min_rec}-{min_abd}/depth_{meandepth}.fa"
    threads: 6
    shell:
        """
        awk 'NR > 1 {{print ">kmer_"NR - 1"\\n"$1}}' {input.kmer_mat} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -num_threads 6 -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -dust no >> {output}
        """

