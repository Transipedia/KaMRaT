#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile f_unsupv-rk-merge.smk --cluster "qsub -q common -l nodes=node27:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr" --jobs 1 -p --latency-wait 60 --rerun-incomplete >> workflow_f_unsuperv-rk-merge_luad.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"

configfile: "luad-config.json"

# Outputs
RES_DIR = config["res_dir"]
MAT_DIR = RES_DIR + "matrices/"
UNSUPV_DIR = RES_DIR + "kamrat_res/all-samples/unsupv-rank-merge/"
IDX_DIR = "/data/work/I2BC/haoliang.xue/kamrat.idx/LUADseo/all-samples/index-norm/"

MTHD_LIST = ["sd", "rsd", "entropy"]

rule all:
    input:
        expand(UNSUPV_DIR + "contig-counts-{mthd}.tsv", mthd = MTHD_LIST)

rule index:
    input:
        MAT_DIR + "kmer-counts.allsmp.tsv"
    output:
        IDX_DIR + "idx-mat.bin",
        IDX_DIR + "idx-meta.bin",
        IDX_DIR + "idx-pos.bin"
    threads: 1
    log:
        IDX_DIR + "log-kamratIndex.allsmp.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat index -intab {input} -outdir {IDX_DIR} -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule unsupv_rank_merge:
    input:
        IDX_DIR + "idx-mat.bin",
        IDX_DIR + "idx-meta.bin",
        IDX_DIR + "idx-pos.bin"
    output:
        rk = UNSUPV_DIR + "kmer-scores-{mthd}.bin",
        mg = UNSUPV_DIR + "contig-counts-{mthd}.tsv"
    params:
        "{mthd}"
    threads: 1
    log:
        UNSUPV_DIR + "log-kamratRank-{mthd}.allsmp.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat rank -idxdir {IDX_DIR} -rankby {params} -outpath {output.rk} -seltop 10000 &> {log}
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {IDX_DIR} -overlap 30-15 -with {output.rk} -interv pearson:0.20 -outpath {output.mg} -withcounts mean &>> {log}
        """

