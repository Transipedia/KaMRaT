#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile a_dataPrepare.smk --cluster "qsub -q common -l nodes=node28:ppn=6" --jobs 3 -p --latency-wait 60 --rerun-incomplete >> workflow_a_dataPrepare.txt &

# Involved programs
CUTADAPT = "/home/haoliang.xue/.conda/envs/cutadapt/bin/cutadapt"
FASTQC = "/home/haoliang.xue/.conda/envs/multiqc/bin/fastqc"
MULTIQC = "/home/haoliang.xue/.conda/envs/multiqc/bin/multiqc"
JELLYFISH = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/jellyfish"
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-da1150c.sif"

#configfile: "luad-config.json"
configfile: "prad-config.json"

# Dataset
FQ_DIR = config["fastq_dir"]
SMP_CONDI_FILE = config["smp_condi_file"]

# Outputs
RES_DIR = config["res_dir"]
TRIM_DIR = RES_DIR + "cutadapt_res/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "matrices/"

# Make sample list from sample-condition file
SMP_LIST = []
with open(SMP_CONDI_FILE) as scf:
    for row in scf:
        SMP_LIST.append(row.split()[0])

# ===== Workflow ===== #
rule all:
    input: 
        expand(MAT_DIR + "sampleshuf.{set}{num}.tsv", set = ["train", "test"], num = range(5)),
        expand(MAT_DIR + "kmer-counts.{set}{num}.tsv", set = ["train", "test"], num = range(5)),
        MAT_DIR + "all-samples.tsv",
        MAT_DIR + "kmer-counts.allsmp.tsv",
        expand(JF_DIR + "{smp}.txt", smp = SMP_LIST),
        TRIM_DIR + "multiqc_report.html"

rule cutadapt:
    input:
        raw1 = FQ_DIR + "{smp}_1.fastq.gz",
        raw2 = FQ_DIR + "{smp}_2.fastq.gz",
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

rule fastqc:
    input:
        trim1 = TRIM_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_DIR + "{smp}_2.trim.fastq.gz"
    output:
        qc1 = TRIM_DIR + "{smp}_1.trim_fastqc.html",
        qc2 = TRIM_DIR + "{smp}_2.trim_fastqc.html"
    log:
        TRIM_DIR + "log-fastqc.{smp}.txt"
    threads: 6
    shell:
        "{FASTQC} --noextract --outdir {TRIM_DIR} --threads 6 {input.trim1} {input.trim2} &> {log}"

rule multiqc:
    input:
        expand(TRIM_DIR + "{smp}_{num}.trim_fastqc.html", smp = SMP_LIST, num = [1, 2])
    output:
        TRIM_DIR + "multiqc_report.html"
    log:
        TRIM_DIR + "log-multiqc.txt"
    threads: 1
    shell:
        "{MULTIQC} --outdir {TRIM_DIR} {TRIM_DIR}"

rule jellyfish:
    input: 
        trim1 = TRIM_DIR + "{smp}_1.trim.fastq.gz",
        trim2 = TRIM_DIR + "{smp}_2.trim.fastq.gz"
    output: 
        jfcount = JF_DIR + "{smp}.jf",
        jfdump = JF_DIR + "{smp}.txt"
    log:
        JF_DIR + "log-jellyfish.{smp}.txt"
    threads: 6
    shell:
        """
        {JELLYFISH} count -m 31 -s 1000000 -C -t 6 -o {output.jfcount} -F 2 <(zcat {input.trim1}) <(zcat {input.trim2}) &>> {log}
        {JELLYFISH} dump -c {output.jfcount} | sort -k 1 -T {JF_DIR} --parallel=6 > {output.jfdump} 2>> {log}
        """

rule split_cv:
    output:
        expand(MAT_DIR + "sampleshuf.{set}{num}.tsv", set = ["train", "test"], num = range(5))
    log:
        MAT_DIR + "log-splitCV.txt"
    threads: 1
    shell:
        """
        export SINGULARITY_BINDPATH="/store:/store"
        singularity exec {KAMRAT_IMG} splitCV.bash -s {SMP_CONDI_FILE} -n 5 -o {MAT_DIR} &> {log}
        """

rule joinCounts:
    input:
        samples = MAT_DIR + "sampleshuf.{set}{num}.tsv",
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        MAT_DIR + "kmer-counts.{set}{num}.tsv"
    log:
        MAT_DIR + "log-joinCounts.{set}{num}.txt"
    threads: 1
    shell:
        """
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {input.samples} > {output} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {input.samples}`
        export SINGULARITY_BINDPATH="/store:/store"
        singularity exec {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {output} 2> {log}
        """

rule joinCounts_all:
    input:
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        dsgn = MAT_DIR + "all-samples.tsv",
        tab = MAT_DIR + "kmer-counts.allsmp.tsv"
    log:
        MAT_DIR + "log-joinCounts.allsmp.txt"
    threads: 1
    shell:
        """
        cp {SMP_CONDI_FILE} {output.dsgn}
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {SMP_CONDI_FILE} > {output.tab} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {SMP_CONDI_FILE}`
        export SINGULARITY_BINDPATH="/store:/store"
        singularity exec {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {output.tab} 2> {log}
        """
