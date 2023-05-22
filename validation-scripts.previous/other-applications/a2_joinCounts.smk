#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile a2_joinCounts.smk --cluster "qsub -q lowprio -l nodes=1:ppn=1 -l mem=200G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_a2_dataPrepare.txt &

# Involved programs
CUTADAPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/cutadapt"
JELLYFISH = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/jellyfish"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Inputs
SMPLIST_FILE = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-list.txt"
SMPCONDI_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-conditions/"
SMP_LIST = []
for f in ["LUADseo.tsv", "PRADtcga.tsv"]:
    with open(SMPCONDI_DIR + f) as scf:
        for row in scf:
            SMP_LIST.append(row.split("\t")[0])

# Outputs
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "joined_matrices/"


# ===== Workflow ===== #
rule all:
    input: 
        expand(MAT_DIR + "{dataset}/sampleshuf.{splitset}{num}.tsv",
               dataset = ["LUADseo-CV", "PRADtcga-CV"], splitset = ["train", "test"], num = range(5)),
        expand(MAT_DIR + "{dataset}/kmer-counts.{splitset}{num}.tsv.gz",
               dataset = ["LUADseo-CV", "PRADtcga-CV"], splitset = ["train", "test"], num = range(5)),
        expand(MAT_DIR + "{x}-smp.tsv", x = ["PRADtcga", "all"]),
        expand(MAT_DIR + "kmer-counts.{x}-smp.tsv.gz", x = ["PRADtcga", "all"])


rule split_cv:
    input:
        SMPCONDI_DIR + "{dataset}.tsv"
    output:
        expand(MAT_DIR + "{ds}-CV/sampleshuf.{splitset}{num}.tsv",
               ds = ["{dataset}"], splitset = ["train", "test"], num = range(5))
    log:
        MAT_DIR + "{dataset}-CV/log-splitCV.txt"
    threads: 1
    params:
        MAT_DIR + "{dataset}-CV/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} splitCV.bash -s {input} -n 5 -o {params} &> {log}
        """

rule joinCounts_cv:
    input:
        samples = MAT_DIR + "{dataset}-CV/sampleshuf.{splitset}{num}.tsv",
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        MAT_DIR + "{dataset}-CV/kmer-counts.{splitset}{num}.tsv.gz"
    log:
        MAT_DIR + "{dataset}-CV/log-joinCounts.{splitset}{num}.txt"
    params:
        MAT_DIR + "{dataset}-CV/kmer-counts.{splitset}{num}.tsv"
    threads: 1
    shell:
        """
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {input.samples} > {params} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {input.samples}`
        apptainer exec -B "/store:/store" {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {params} 2> {log}
        gzip {params}
        """

rule joinCounts_datasets:
    input:
        smpcondi_file = SMPCONDI_DIR + "{dataset}.tsv",
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        dsgn = MAT_DIR + "{dataset}-smp.tsv",
        tab = MAT_DIR + "kmer-counts.{dataset}-smp.tsv.gz"
    log:
        MAT_DIR + "log-joinCounts.{dataset}-smp.txt"
    params:
        MAT_DIR + "kmer-counts.{dataset}-smp.tsv"
    threads: 1
    shell:
        """
        cp {input.smpcondi_file} {output.dsgn}
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {output.dsgn} > {params} # print the header row 
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {output.dsgn}` 
        apptainer exec -B "/store:/store" {KAMRAT_IMG} joinCounts -r 1 -a 1 ${{path_list[@]}} >> {params} 2> {log}
        gzip {params}
        """

rule joinCounts_allsmp:
    input:
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        dsgn = MAT_DIR + "all-smp.tsv",
        tab = MAT_DIR + "kmer-counts.all-smp.tsv.gz"
    log:
        MAT_DIR + "log-joinCounts.all-smp.txt"
    params:
        tab = MAT_DIR + "kmer-counts.all-smp.tsv"
    threads: 1
    shell:
        """
        cat {SMPCONDI_DIR}[LP]*.tsv > {output.dsgn}
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {output.dsgn} > {params} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {output.dsgn}`
        apptainer exec -B "/store:/store" {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {params} 2> {log}
        gzip {params}
        """
