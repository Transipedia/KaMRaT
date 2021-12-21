#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile b_kamrat.smk --cluster "qsub -q common -l nodes=node13:ppn=1" --jobs 2 -p --latency-wait 60 --rerun-incomplete >> workflow_b_kamrat.txt &

# Involved programs
RSCRIPT = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"
BLASTN = "/usr/bin/blastn"

# Inputs
REF_FASTA = "/store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa"
BLAST_DB = "/store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa"

# Outputs
RES_DIR = "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-merge/"
MAT_DIR = RES_DIR + "matrices/"
IDX_DIR = "/data/work/I2BC/haoliang.xue/kamrat.idx/polyester_simu/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"
BLASTN_DIR = RES_DIR + "blastn_res/"

P_LIST = [i * 10 for i in range(4, 11)]
INTERV_LIST = ["pearson", "spearman", "mac"]
THRES_LIST = [str(i / 10) for i in range(1, 11)]

# ===== Workflow ===== #
rule all:
    input:
        expand(KAMRAT_DIR + "{p}pct/merged-counts.{interv}_{thres}.tsv", p = P_LIST, interv = INTERV_LIST, thres = THRES_LIST),
        expand(KAMRAT_DIR + "{p}pct/merged-counts.none.tsv", p = P_LIST),
        expand(BLASTN_DIR + "{p}pct/ctg-aligned.none.tsv", p = P_LIST),
        expand(BLASTN_DIR + "{p}pct/ctg-aligned.{interv}_{thres}.tsv", p = P_LIST, interv = INTERV_LIST, thres = THRES_LIST)

rule kamrat_index:
    input:
        MAT_DIR + "kmer-counts-{p}pct.tsv"
    output:
        IDX_DIR + "index-{p}pct/idx-meta.bin",
        IDX_DIR + "index-{p}pct/idx-pos.bin",
        IDX_DIR + "index-{p}pct/idx-mat.bin"
    log:
        IDX_DIR + "log-kamrat_index-{p}pct.txt"
    threads: 1
    params:
        "{p}"
    shell:
        """
        export SINGULARITY_BINDPATH="/store:/store,/data:/data"
        mkdir -p {IDX_DIR}index-{params}pct/
        singularity exec {KAMRAT_IMG} kamrat index -intab {input} -outdir {IDX_DIR}index-{params}pct -klen 31 -unstrand -nfbase 1000000000 &> {log}
        """

rule kamrat_merge_none:
    input:
        IDX_DIR + "index-{p}pct/idx-meta.bin",
        IDX_DIR + "index-{p}pct/idx-pos.bin",
        IDX_DIR + "index-{p}pct/idx-mat.bin"
    output:
        KAMRAT_DIR + "{p}pct/merged-counts.none.tsv"
    threads: 1
    params:
        "{p}"
    log:
        KAMRAT_DIR + "{p}pct/log-kamrat-merge-none.txt"
    shell:
        """
        export SINGULARITY_BINDPATH="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {IDX_DIR}index-{params}pct -overlap 30-15 -interv none -outpath {output} -withcounts rep &> {log}
        """

rule kamrat_merge_interv:
    input:
        IDX_DIR + "index-{p}pct/idx-meta.bin",
        IDX_DIR + "index-{p}pct/idx-pos.bin",
        IDX_DIR + "index-{p}pct/idx-mat.bin"
    output:
        KAMRAT_DIR + "{p}pct/merged-counts.{interv}_{thres}.tsv"
    log:
        KAMRAT_DIR + "{p}pct/log-kamrat-merge-{interv}_{thres}.txt"
    threads: 1
    params:
        p = "{p}",
        interv = "{interv}",
        thres = "{thres}"
    shell:
        """
        export SINGULARITY_BINDPATH="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {IDX_DIR}index-{params.p}pct -overlap 30-15 -interv {params.interv}:{params.thres} \
                                                   -outpath {output} -withcounts rep &> {log}
        """

rule blastn:
    input:
        KAMRAT_DIR + "{p}pct/merged-counts.{mode}.tsv"
    output:
        BLASTN_DIR + "{p}pct/ctg-aligned.{mode}.tsv"
    threads: 1
    params:
        BLASTN_DIR + "{p}pct/ctg-seq.{mode}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {BLAST_DB} -query {params} -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 1 -dust no >> {output}
        """
