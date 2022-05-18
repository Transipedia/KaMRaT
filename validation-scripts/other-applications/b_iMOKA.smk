#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile b_iMOKA.smk --cluster "qsub -l nodes=node29:ppn=5 -q ssfa -m ea -M haoliang.xue@i2bc.paris-saclay.fr" --jobs 1 -p --latency-wait 60 --rerun-incomplete >> workflow_b_iMOKA_luad.txt &

# Involved programs
IMOKA_IMG = "/home/haoliang.xue/tools/iMOKA-1.1.img"

# configfile: "luad-config.json"
configfile: "prad-config.json"

# Outputs
RES_DIR = config["res_dir"]
JF_DIR = RES_DIR + "jellyfish_res/"
SMPINFO_DIR = RES_DIR + "matrices/"
IMOKA_DIR = RES_DIR + "imoka_res_check/"


rule all:
    input:
        expand(IMOKA_DIR + "CV{num}/{mode}/prediction_test.txt", num = range(5), mode = ["simple", "default"])


rule imoka_create_all: # first, prepare all sample information for one-turn index across all samples
    input:
        train0 = SMPINFO_DIR + "sampleshuf.train0.tsv",
        test0 = SMPINFO_DIR + "sampleshuf.test0.tsv"
    output:
        imoka_in = IMOKA_DIR + "input_allsmp",
        imoka_json = IMOKA_DIR + "matrix_train.json"
    threads: 5
    log:
        IMOKA_DIR + "log-imoka-createall.txt"
    shell:
        """
        awk -v jfdir="{JF_DIR}" '$1 != "sample" {{print jfdir"/"$1".txt\t"$1"\t"$2}}' {input.train0} {input.test0} > {output.imoka_in}

        export SINGULARITY_BINDPATH="/store:/store"
        export OMP_NUM_THREADS=5
        singularity exec {IMOKA_IMG} iMOKA_core create -i {output.imoka_in} -o {output.imoka_json} &> {log}
        """

rule imoka_create_byCV: # real iMOKA create, fold by fold
    input:
        smpinfo = SMPINFO_DIR + "sampleshuf.{set}{num}.tsv",
        tmp = IMOKA_DIR + "matrix_train.json"
    output:
        imoka_in = IMOKA_DIR + "CV{num}/input_{set}",
        imoka_json = IMOKA_DIR + "CV{num}/matrix_{set}.json"
    threads: 5
    log:
        IMOKA_DIR + "log-imoka-create-{set}{num}.txt"
    shell:
        """
        awk -v jfdir="{JF_DIR}" '$1 != "sample" {{print jfdir"/"$1".txt\t"$1"\t"$2}}' {input.smpinfo} > {output.imoka_in}

        export SINGULARITY_BINDPATH="/store:/store"
        export OMP_NUM_THREADS=5
        singularity exec {IMOKA_IMG} iMOKA_core create -i {output.imoka_in} -o {output.imoka_json} &> {log}
        """

rule imoka_reduce_aggregate_norep: # by CV, no repetition
    input:
        IMOKA_DIR + "CV{num}/matrix_train.json"
    output:
        red = IMOKA_DIR + "CV{num}/simple/reduced.matrix",
        aggmat = IMOKA_DIR + "CV{num}/simple/aggregated.kmers.matrix"
    threads: 5
    params:
        agg = IMOKA_DIR + "CV{num}/simple/",
    resources:
        mem_mb = 120 * 1024 # 120G memory
    log:
        IMOKA_DIR + "CV{num}/simple/log-imoka-red-agg.txt"
    shell:
        """
        export IMOKA_MAX_MEM_GB=100
        export SINGULARITY_BINDPATH="/store:/store"
        export OMP_NUM_THREADS=5
        singularity exec {IMOKA_IMG} iMOKA_core reduce -c 1 -i {input} -o {output.red} &> {log}
        singularity exec {IMOKA_IMG} iMOKA_core aggregate -i {output.red} -o {params.agg}aggregated -c {input} -m nomap &>> {log}
        """

rule imoka_reduce_aggregate_default: # by CV, default 100 repetition
    input:
        IMOKA_DIR + "CV{num}/matrix_train.json"
    output:
        red = IMOKA_DIR + "CV{num}/default/reduced.matrix",
        aggmat = IMOKA_DIR + "CV{num}/default/aggregated.kmers.matrix"
    threads: 5
    params:
        agg = IMOKA_DIR + "CV{num}/default/",
    resources:
        mem_mb = 120 * 1024 # 120G memory
    log:
        IMOKA_DIR + "CV{num}/default/log-imoka-red-agg.txt"
    shell:
        """
        export IMOKA_MAX_MEM_GB=100
        export SINGULARITY_BINDPATH="/store:/store"
        export OMP_NUM_THREADS=5
        singularity exec {IMOKA_IMG} iMOKA_core reduce -i {input} -o {output.red} &> {log}
        singularity exec {IMOKA_IMG} iMOKA_core aggregate -i {output.red} -o {params.agg}aggregated -c {input} -m nomap &>> {log}
        """

rule build_and_test: # by CV
    input:
        agg = IMOKA_DIR + "CV{num}/{mode}/aggregated.kmers.matrix",
        test_json = IMOKA_DIR + "CV{num}/matrix_test.json"
    output:
        rf = directory(IMOKA_DIR + "CV{num}/{mode}/randomforest_models"),
        test_mat = IMOKA_DIR + "CV{num}/{mode}/extracted-test.mat",
        pred = IMOKA_DIR + "CV{num}/{mode}/prediction_test.txt"
    threads: 5
    params:
        rf_out = IMOKA_DIR + "CV{num}/{mode}/randomforest"
    log:
        IMOKA_DIR + "CV{num}/{mode}/log-imoka-build-and-test.txt"
    shell:
        """
        export IMOKA_MAX_MEM_GB=100
        export SINGULARITY_BINDPATH="/store:/store"
        export OMP_NUM_THREADS=5
        singularity exec {IMOKA_IMG} random_forest.py -m 100 {input.agg} {params.rf_out} &> {log}
        singularity exec {IMOKA_IMG} iMOKA_core extract -i {output.rf}/0_RF.pickle.features -o {output.test_mat} -s {input.test_json} &>> {log}
        singularity exec {IMOKA_IMG} predict.py {output.test_mat} {output.rf}/0_RF.pickle {output.pred} &>> {log}
        """
