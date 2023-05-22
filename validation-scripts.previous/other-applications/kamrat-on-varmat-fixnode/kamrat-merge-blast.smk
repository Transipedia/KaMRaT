#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile kamrat-merge-blast.smk --cluster "qsub -q lowprio -l nodes=1:ppn=1 -l mem=200G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_kamrat-merge-blast.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"

# Inputs and outputs
BLAST_DB = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/blast_db/"
REF_PATH = BLAST_DB + "gencode.v34.transcripts.fa"
KAMRAT_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/kamrat_on_varmat/apply-on-varmat/140-smp/"

rule all:
    input:
        expand(KAMRAT_DIR + "contig-alignment.{mode}.tsv", mode = ["mac", "pearson", "spearman", "none"])

rule makeblastdb:
    input:
        REF_PATH
    output:
        BLAST_DB + "gencode.v34.transcripts.fa",
        BLAST_DB + "gencode.v34.transcripts.fa.nhr",
        BLAST_DB + "gencode.v34.transcripts.fa.nin",
        BLAST_DB + "gencode.v34.transcripts.fa.nog",
        BLAST_DB + "gencode.v34.transcripts.fa.nsd",
        BLAST_DB + "gencode.v34.transcripts.fa.nsi",
        BLAST_DB + "gencode.v34.transcripts.fa.nsq"
    threads: 1
    shell:
        """
        {MKBLASTDB} -in {input} -dbtype nucl -parse_seqids -out {BLAST_DB}gencode.v34.transcripts.fa
        """

rule blastn:
    input:
        ctg_file = KAMRAT_DIR + "kamrat-index-merge.{mode}.tsv",
        blast_db = REF_PATH
    output:
        KAMRAT_DIR + "contig-alignment.{mode}.tsv"
    threads: 1
    params:
        KAMRAT_DIR + "contig-seq.{mode}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input.ctg_file} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid\\tqcovs\\tqacc" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -max_hsps 1 -max_target_seqs 1 -evalue 0.001 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid qcovs qacc" -num_threads 1 -dust no >> {output}
        """
