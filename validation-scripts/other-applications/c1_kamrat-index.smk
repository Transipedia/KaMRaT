#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile c1_kamrat-index.smk --cluster "qsub -q lowprio -l nodes=1:ppn=1 -l mem=200G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_c1_kamrat-index.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Inputs and outputs
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "joined_matrices/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"
DATASET = "LUADseo"  # LUADseo, PRADtcga, 20smp, 50smp, 80smp, 110smp, 140smp, 170smp, 200smp, 230smp

rule all:
    input:
        expand(KAMRAT_DIR + "{dataset}/CV{num}/index-{cvset}/{f}.bin",
               dataset = ["LUADseo", "PRADtcga"], num = range(5),
               cvset = ["train", "test"], f = ["idx-meta", "idx-pos", "idx-mat"]),
        expand(KAMRAT_DIR + "{dataset}/allsmp/index/{f}.bin",
               dataset = ["LUADseo", "PRADtcga"], f = ["idx-meta", "idx-pos", "idx-mat"]),
        expand(KAMRAT_DIR + "{num}smp/index/{f}.bin",
               num = [20, 50, 80, 110, 140, 170, 200, 230], f = ["idx-meta", "idx-pos", "idx-mat"])


rule kamrat_index_dataset_cv:
    input:
        MAT_DIR + "{dataset}/kmer-counts.{cvset}{num}.tsv.gz"
    output:
        expand(KAMRAT_DIR + "{dataset}/CV{num}/index-{cvset}/{f}.bin",
               dataset = ["{dataset}"], num = ["{num}"], cvset = ["{cvset}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "{dataset}/CV{num}/index-{cvset}/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "{dataset}/CV{num}/index-{cvset}/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                    -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_index_dataset_all:
    input:
        MAT_DIR + "{dataset}/kmer-counts.allsmp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "{dataset}/allsmp/index/{f}.bin",
               dataset = ["{dataset}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "{dataset}/allsmp/index/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "{dataset}/allsmp/index/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                    -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_index_varsize:
    input:
        MAT_DIR + "kmer-counts_varsize.{num}smp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "{num}smp/index/{f}.bin",
               num = ["{num}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "{num}smp/index/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "{num}smp/index/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                    -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """
