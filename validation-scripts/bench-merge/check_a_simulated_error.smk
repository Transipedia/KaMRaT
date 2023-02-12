#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile check_a_simulated_error.smk --cluster "qsub -q common -l nodes=1:ppn=6 -l mem=20G -l walltime=300:00:00" --jobs 8 -p --latency-wait 60 --rerun-incomplete >> workflow_check_a_simulated_error.txt &


# Parameters
SMP_LIST = ["sample_01", "sample_02", "sample_03", "sample_04", "sample_05", "sample_06", "sample_07", "sample_08", "sample_09", "sample_10",
            "sample_11", "sample_12", "sample_13", "sample_14", "sample_15", "sample_16", "sample_17", "sample_18", "sample_19", "sample_20"]
MEAN_DEPTH_LIST = list(range(1, 11))

# Involved programs
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"

# Inputs and outputs
REF_FASTA = "/store/plateformes/CALCUL/SSFA_KaMRaT/Data/gc34.50-53.fa"
RES_DIR = f"/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
BLAST_DB = RES_DIR + "blast_db/"
PLSTR_DIR = RES_DIR + "polyester_res/"
CHECK_READ_DIR = RES_DIR + "check_read_error/"


# ===== Workflow ===== #
rule all:
    input:
        expand(CHECK_READ_DIR + "{mode}/depth_{meandepth}/{smp}_{num}.align.tsv",
               mode = ["err-free", "err-illumina5"], smp = SMP_LIST, num = [1, 2], meandepth = MEAN_DEPTH_LIST)

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
        read_file = PLSTR_DIR + "{mode}/depth_{meandepth}/{file}.fasta",
        blast_db = BLAST_DB + "gc34.50-53.fa"
    output:
        CHECK_READ_DIR + "{mode}/depth_{meandepth}/{file}.align.tsv"
    threads: 6
    shell:
        """
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {input.blast_db} -query {input.read_file} -num_threads 6 -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -dust no >> {output}
        """

