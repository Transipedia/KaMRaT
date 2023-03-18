#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile b2_kamrat-on-top-kmers.smk --cluster "qsub -q common -l nodes=1 -l mem=100G -l walltime=5:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 5 -p --latency-wait 60 --rerun-incomplete >> workflow_b2_kamrat-on-top.txt &

# Involved programs
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"

# Inputs
DSGN_FILE = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/polyester_res/depth_10/sim_rep_info.txt"
RES_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
BLASTDB_DIR = RES_DIR + "blast_db/"

# Outputs
MAT_DIR = RES_DIR + "matrices/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"

INTERV_LIST = ["pearson", "spearman", "mac"]
THRES_LIST = [str(i / 10) for i in range(1, 11)]
TOP_PCT_LIST = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# ===== Workflow ===== #
rule all:
    input:
        expand(KAMRAT_DIR + "depth_10/ctg-aligned.none.topkmers-ttest.padj-{top_pct}.tsv", 
               top_pct = TOP_PCT_LIST),
        expand(KAMRAT_DIR + "depth_10/ctg-aligned.{interv}_0.2.topkmers-ttest.padj-{top_pct}.tsv",
               interv=INTERV_LIST, top_pct = TOP_PCT_LIST),
        expand(KAMRAT_DIR + "depth_10/ctg-aligned.{interv}_{thres}.topkmers-ttest.padj-0.4.tsv",
               interv=INTERV_LIST, thres=THRES_LIST)


rule kamrat_index:
    input:
        MAT_DIR + "depth_10/kmer-counts.tsv"
    output:
        KAMRAT_DIR + "depth_10/index/idx-meta.bin",
        KAMRAT_DIR + "depth_10/index/idx-pos.bin",
        KAMRAT_DIR + "depth_10/index/idx-mat.bin"
    params:
        KAMRAT_DIR + "depth_10/index/"
    log:
        KAMRAT_DIR + "depth_10/index/log-kamrat_index.txt"
    threads: 1
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat index -intab {input} \
                                                                  -outdir {params} \
                                                                  -klen 31 -unstrand -nfbase 1000000000 &> {log}
        """

rule kamrat_rank:
    input:
        design = KAMRAT_DIR + "design.txt",
        idx1 = KAMRAT_DIR + "depth_10/index/idx-meta.bin",
        idx2 = KAMRAT_DIR + "depth_10/index/idx-pos.bin",
        idx3 = KAMRAT_DIR + "depth_10/index/idx-mat.bin"
    output:
        KAMRAT_DIR + "depth_10/top-kmers-{top_pct}.ttest.padj.bin"
    params:
        idxdir = KAMRAT_DIR + "depth_10/index/",
        top_pct = "{top_pct}"
    threads: 1
    log:
        KAMRAT_DIR + "depth_10/log-kamrat-rank-{top_pct}.txt"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat rank -idxdir {params.idxdir} \
                                                                 -rankby ttest.padj \
                                                                 -design {input.design} \
                                                                 -outpath {output} \
                                                                 -seltop {params.top_pct} &> {log}
        """

rule kamrat_merge_none:
    input:
        topres = KAMRAT_DIR + "depth_10/top-kmers-{top_pct}.ttest.padj.bin",
        idx1 = KAMRAT_DIR + "depth_10/index/idx-meta.bin",
        idx2 = KAMRAT_DIR + "depth_10/index/idx-pos.bin",
        idx3 = KAMRAT_DIR + "depth_10/index/idx-mat.bin"
    output:
        KAMRAT_DIR + "depth_10/merged-counts.none.topkmers-ttest.padj-{top_pct}.tsv"
    params:
        KAMRAT_DIR + "depth_10/index/",
    threads: 1
    log:
        KAMRAT_DIR + "depth_10/log-kamrat-merge-none.topkmers-ttest.padj-{top_pct}.txt"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params} \
                                                                  -overlap 30-15 -interv none \
                                                                  -with {input.topres} \
                                                                  -outpath {output} -withcounts rep &> {log}
        """

rule kamrat_merge_interv:
    input:
        topres = KAMRAT_DIR + "depth_10/top-kmers-{top_pct}.ttest.padj.bin",
        idx1 = KAMRAT_DIR + "depth_10/index/idx-meta.bin",
        idx2 = KAMRAT_DIR + "depth_10/index/idx-pos.bin",
        idx3 = KAMRAT_DIR + "depth_10/index/idx-mat.bin"
    output:
        KAMRAT_DIR + "depth_10/merged-counts.{interv}_{thres}.topkmers-ttest.padj-{top_pct}.tsv"
    log:
        KAMRAT_DIR + "depth_10/log-kamrat-merge-{interv}_{thres}.topkmers-ttest.padj-{top_pct}.txt"
    threads: 1
    params:
        idxdir = KAMRAT_DIR + "depth_10/index/",
        interv = "{interv}",
        thres = "{thres}"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params.idxdir} \
                                                                  -overlap 30-15 -interv {params.interv}:{params.thres} \
                                                                  -with {input.topres} \
                                                                  -outpath {output} -withcounts rep &> {log}
        """

rule blastn:
    input:
        ctg_file = KAMRAT_DIR + "depth_10/merged-counts.{mode}.topkmers-ttest.padj-{top_pct}.tsv",
        blast_db = BLASTDB_DIR + "gc34.50-53.fa"
    output:
        KAMRAT_DIR + "depth_10/ctg-aligned.{mode}.topkmers-ttest.padj-{top_pct}.tsv"
    threads: 1
    params:
        KAMRAT_DIR + "depth_10/ctg-seq.{mode}.topkmers-ttest.padj-{top_pct}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input.ctg_file} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 1 -dust no >> {output}
        """
