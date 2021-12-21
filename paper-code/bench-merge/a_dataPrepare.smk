#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile a_dataPrepare.smk --cluster "qsub -q common -l nodes=node22:ppn=6" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_a_dataPrepare.txt &

# Involved programs
POLYESTER = "/store/USERS/haoliang.xue/development/KaMRaT/paper-code/bench-merge/polyester.R"
RSCRIPT = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript"
JELLYFISH = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/jellyfish"
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"

# Inputs
REF_FASTA = "/store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa"

# Outputs
RES_DIR = "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-merge/"
PLSTR_DIR = RES_DIR + "polyester_res/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "matrices/"

SMP_LIST = ["sample_01", "sample_02", "sample_03", "sample_04", "sample_05", "sample_06", "sample_07", "sample_08", "sample_09", "sample_10",
            "sample_11", "sample_12", "sample_13", "sample_14", "sample_15", "sample_16", "sample_17", "sample_18", "sample_19", "sample_20"]

# ===== Workflow ===== #
rule all:
    input:
        expand(MAT_DIR + "kmer-counts-{p}pct.tsv", p = [i * 10 for i in range(4, 11)])

rule polyester:
    output:
        expand(PLSTR_DIR + "{smp}_{num}.fasta", smp = SMP_LIST, num = [1, 2])
    log:
        PLSTR_DIR + "log-polyester.txt"
    shell:
        """
        {RSCRIPT} {POLYESTER} {REF_FASTA} {PLSTR_DIR} 2> {log}
        """

rule jellyfish:
    input:
        fa1 = PLSTR_DIR + "{smp}_1.fasta",
        fa2 = PLSTR_DIR + "{smp}_2.fasta"
    output:
        JF_DIR + "{smp}.txt"
    params:
        JF_DIR + "{smp}.jf"
    log:
        JF_DIR + "log-jellyfish-{smp}.txt"
    shell:
        """
        {JELLYFISH} count -m 31 -s 1000000 -t 6 -o {params} -F 2 -C {input.fa1} {input.fa2} &> {log}
	    {JELLYFISH} dump -c {params} | sort -k 1 -T {JF_DIR} --parallel=6 > {output} 2>> {log}
        rm -f {params}
        """

rule joinCounts:
    input:
        expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        MAT_DIR + "kmer-counts-100pct.tsv"
    log:
        MAT_DIR + "log-joinCounts.txt"
    shell:
        """
        echo -en "tag\\t" > {output}
        echo `ls {input} | sed 's/.txt//g' | sed 's|.*/sample|sample|g'` | sed 's/ /\\t/g' >> {output}
        export SINGULARITY_BINDPATH="/store:/store"
        singularity exec {KAMRAT_IMG} joinCounts -r 1 -a 1 {input} >> {output} 2> {log}
        """

rule gene_subtab:
    input:
        MAT_DIR + "kmer-counts-100pct.tsv"
    output:
        MAT_DIR + "kmer-counts-{p}pct.tsv"
    params:
        "{p}"
    shell:
        """
        head -n1 {input} > {output}
        nline=$[`wc -l {input} | awk '{{print $1}}'` - 1]
        sed '1d' {input} | shuf -n $[$nline * {params} / 100] >> {output}
        """
