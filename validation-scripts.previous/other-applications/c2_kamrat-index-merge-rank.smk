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

rule all:
    input:
        expand(KAMRAT_DIR + "{dataset}/merged-kmer-counts.{mode}.tsv",
               dataset = ["LUADseo/allsmp", "PRADtcga/allsmp"],
               mode = ["mac", "pearson", "spearman", "none"])

rule split_mat:
    input:
        tab = MAT_DIR + "kmer-counts.all-smp.tsv.gz",
        meta = MAT_DIR + "all-smp.tsv"
    output:
        tab = MAT_DIR + "kmer-counts.{n}-smp.tsv.gz",
        meta = MAT_DIR + "{n}-smp.tsv"
    log:
        MAT_DIR + "log-splitmat.{n}-smp.txt"
    params:
        nsmp = "{n}",
        tab = MAT_DIR + "kmer-counts.{n}-smp.tsv"
    threads: 1
    shell:
        """
        head -n$(({params.nsmp} + 1)) {input.meta} > {output.meta}
        zcat {input.tab} | awk -v nsmp="{params.nsmp}" '{{for (i = 1; i < nsmp; i++) printf("%s\\t", $i); print $nsmp}}' > {params.tab} 2> {log}
        gzip {params.tab}
        """

rule kamrat_index:
    input:
        MAT_DIR + "kmer-counts.{x}-smp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "{x}-smp/index/{f}.bin", x = ["{x}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "{x}-smp/index/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "{x}-smp/index/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                    -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_merge:
    input:
        KAMRAT_DIR + "{x}-smp/index/"
    output:
        KAMRAT_DIR + "{x}-smp/merged-contigs.tsv"
    threads: 1
    resources:
        mem_mb = 250 * 1024 # 120G memory
    log:
        KAMRAT_DIR + "{x}-smp/log-kamrat-merge.txt"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat merge -idxdir {input} -overlap 30-15 \
                                                                    -interv pearson:0.20 -outpath {output} &> {log}
        """
