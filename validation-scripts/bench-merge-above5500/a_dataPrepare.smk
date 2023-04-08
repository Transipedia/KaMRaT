#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile a_dataPrepare.smk --cluster "qsub -q common -l nodes=node13:ppn=6 -l mem=30G -l walltime=300:00:00" --jobs 10 -p --latency-wait 60 --rerun-incomplete >> workflow_a_dataPrepare.txt &

# Involved programs
MUTATIONSIMULATOR = "/home/haoliang.xue_ext/miniconda3/envs/py310/bin/mutation-simulator"
POLYESTER = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/bench-merge-new/polyester.R"
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
JELLYFISH = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/jellyfish"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Inputs
REF_FASTA = "/store/plateformes/CALCUL/SSFA_KaMRaT/Data/gc34.above5500nc.fa"

# Parameters
MEAN_DEPTH_LIST = [0.05, 0.2, 0.4, 0.6, 0.8, 1]

# Outputs
RES_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge-above5500/"
SMPREF_DIR = RES_DIR + "smpref_res/"
PLSTR_DIR = RES_DIR + "polyester_res/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "matrices/"

# ===== Workflow ===== #
rule all:
    input:
        expand(MAT_DIR + "err-free/depth_{meandepth}/kmer-counts-1-1.tsv",
               meandepth = MEAN_DEPTH_LIST),

rule simulate_mutations:
    input:
        REF_FASTA
    output:
        SMPREF_DIR + "ref{iref}_ms.fa"
    params:
        iref = "{iref}",
        out = SMPREF_DIR + "ref{iref}"
    log:
        SMPREF_DIR + "log-mutation-simulator_{iref}.txt"
    shell:
        """
        if [ {params.iref} -le 5 ]
        then
            {MUTATIONSIMULATOR} -o {params.out} {input} args -sn 0.001 -in 0.001 -de 0.001 -du 0.001 -iv 0.001 -tl 0.001 -inmax 10 -demax 10 -dumax 10 -ivmax 10 -tlmax 10
        else
            {MUTATIONSIMULATOR} -o {params.out} {input} args -sn 0.005 -in 0.005 -de 0.005 -du 0.005 -iv 0.005 -tl 0.005 -inmax 10 -demax 10 -dumax 10 -ivmax 10 -tlmax 10
        fi
        echo "" >> {output}
        """

rule polyester:
    input:
        SMPREF_DIR + "ref{iref}_ms.fa"
    output:
        expand(PLSTR_DIR + "{mode}/depth_{meandepth}/ref{iref}_ms/{smp}_{num}.fasta",
               mode = ["{mode}"], meandepth = ["{meandepth}"], iref = ["{iref}"],
               smp = ["sample_01", "sample_02"], num = [1, 2])
    params:
        meandepth = "{meandepth}",
        mode = "{mode}"
    log:
        PLSTR_DIR + "{mode}/depth_{meandepth}/ref{iref}_ms/log-polyester.txt"
    shell:
        """
        {RSCRIPT} {POLYESTER} {input} {PLSTR_DIR} {params.meandepth} {params.mode} 2> {log}
        """

rule jellyfish:
    input:
        fa1 = PLSTR_DIR + "{mode}/depth_{meandepth}/ref{iref}_ms/{smp}_1.fasta",
        fa2 = PLSTR_DIR + "{mode}/depth_{meandepth}/ref{iref}_ms/{smp}_2.fasta"
    output:
        JF_DIR + "{mode}/depth_{meandepth}/ref{iref}-{smp}.txt"
    params:
        JF_DIR + "{mode}/depth_{meandepth}/ref{iref}-{smp}.jf"
    log:
        JF_DIR + "{mode}/depth_{meandepth}/log-jellyfish-ref{iref}-{smp}.txt"
    shell:
        """
        {JELLYFISH} count -m 31 -s 1000000 -t 6 -o {params} -F 2 -C {input.fa1} {input.fa2} &> {log}
	    {JELLYFISH} dump -c {params} | sort -k 1 > {output} 2>> {log}
        rm -f {params}
        """

rule joinCounts:
    input:
        expand(JF_DIR + "{mode}/depth_{meandepth}/ref{iref}-{smp}.txt",
               mode = ["{mode}"], meandepth = ["{meandepth}"],
               iref = list(range(1, 11)), smp = ["sample_01", "sample_02"])
    output:
        MAT_DIR + "{mode}/depth_{meandepth}/kmer-counts-{min_rec}-{min_abd}.tsv"
    params:
        min_rec = "{min_rec}",
        min_abd = "{min_abd}"
    log:
        MAT_DIR + "{mode}/depth_{meandepth}/log-joinCounts.txt"
    shell:
        """
        echo -en "tag\\t" > {output}
        echo `ls {input} | sed 's/.txt//g' | sed 's|.*/sample|sample|g'` | sed 's/ /\\t/g' >> {output}
        apptainer exec -B "/data:/data" {KAMRAT_IMG} joinCounts -r {params.min_rec} -a {params.min_abd} {input} >> {output} 2> {log}
        """
