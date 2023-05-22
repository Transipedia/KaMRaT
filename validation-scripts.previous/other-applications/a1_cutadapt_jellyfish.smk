#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile a1_cutadapt_jellyfish.smk --cluster "qsub -q common -l nodes=1:ppn=6 -l mem=30G -l walltime=30:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_a_dataPrepare.txt &

# Involved programs
CUTADAPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/cutadapt"
JELLYFISH = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/jellyfish"

# Inputs
SMP_LIST_FILE = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-list.txt"
SMP_LIST = []
with open(SMP_LIST_FILE) as scf:
    for row in scf:
        SMP_LIST.append(row.split(".")[0][:-2])
SMP_LIST = list(set(SMP_LIST))

# Outputs
FASTQ_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Data/real_fastq/"
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
TRIM_DIR = RES_DIR + "trimmed_res/" 
JF_DIR = RES_DIR + "jellyfish_res/"

# ===== Workflow ===== #
rule all:
    input: 
        expand(JF_DIR + "{smp}.txt", smp = SMP_LIST),

rule cutadapt:
    input:
        raw1 = FASTQ_DIR + "{smp}_1.fastq.gz",
        raw2 = FASTQ_DIR + "{smp}_2.fastq.gz",
    output:
        trim1 = TRIM_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_DIR + "{smp}_2.trim.fastq.gz"
    log:
        TRIM_DIR + "log-cutadapt.{smp}.txt"
    threads: 6
    shell:
        """
        {CUTADAPT} -j 6 -q 12,12 -m 31 -o {output.trim1} -p {output.trim2} \
                   -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                   {input.raw1} {input.raw2} &> {log}
        """

rule jellyfish:
    input: 
        trim1 = TRIM_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_DIR + "{smp}_2.trim.fastq.gz"
    output: 
        jfdump = JF_DIR + "{smp}.txt"
    log:
        JF_DIR + "log-jellyfish.{smp}.txt"
    params:
        JF_DIR + "{smp}.jf"
    threads: 6
    shell:
        """
        {JELLYFISH} count -m 31 -s 1000000 -C -t 6 -o {params} -F 2 <(zcat {input.trim1}) <(zcat {input.trim2}) &>> {log}
        {JELLYFISH} dump -c {params} | sort -k 1 > {output.jfdump} 2>> {log}
        rm -f {params}
        """
