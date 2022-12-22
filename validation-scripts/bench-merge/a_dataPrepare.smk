#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile a_dataPrepare.smk --cluster "qsub -q common -l nodes=node15 -l mem=20G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 3 -p --latency-wait 60 --rerun-incomplete >> workflow_a_dataPrepare.txt &

# Involved programs
POLYESTER = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/bench-merge/polyester.R"
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
JELLYFISH = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/jellyfish"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Inputs
REF_FASTA = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Data/bench-merge/gencode.v34.transcripts.fa"

# Outputs
RES_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
PLSTR_DIR = RES_DIR + "polyester_res/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "matrices/"

SMP_LIST = ["sample_01", "sample_02", "sample_03", "sample_04", "sample_05", "sample_06", "sample_07", "sample_08", "sample_09", "sample_10",
            "sample_11", "sample_12", "sample_13", "sample_14", "sample_15", "sample_16", "sample_17", "sample_18", "sample_19", "sample_20"]
MEAN_DEPTH_LIST = [0.1, 0.2, 0.5, 1, 2, 5, 10]

# ===== Workflow ===== #
rule all:
    input:
        expand(MAT_DIR + "depth_{meandepth}/kmer-counts.tsv", meandepth = MEAN_DEPTH_LIST)

rule polyester:
    output:
        expand(PLSTR_DIR + "depth_{meandepth}/{smp}_{num}.fasta", smp = SMP_LIST, num = [1, 2], meandepth = ["{meandepth}"])
    params:
        meandepth = "{meandepth}"
    log:
        PLSTR_DIR + "depth_{meandepth}/log-polyester.txt"
    shell:
        """
        {RSCRIPT} {POLYESTER} {REF_FASTA} {PLSTR_DIR} {params.meandepth} 2> {log}
        """

rule jellyfish:
    input:
        fa1 = PLSTR_DIR + "depth_{meandepth}/{smp}_1.fasta",
        fa2 = PLSTR_DIR + "depth_{meandepth}/{smp}_2.fasta"
    output:
        JF_DIR + "depth_{meandepth}/{smp}.txt"
    params:
        JF_DIR + "depth_{meandepth}/{smp}.jf"
    log:
        JF_DIR + "depth_{meandepth}/log-jellyfish-{smp}.txt"
    shell:
        """
        {JELLYFISH} count -m 31 -s 1000000 -t 3 -o {params} -F 2 -C {input.fa1} {input.fa2} &> {log}
	    {JELLYFISH} dump -c {params} | sort -k 1 -T {JF_DIR} --parallel=3 > {output} 2>> {log}
        rm -f {params}
        """

rule joinCounts:
    input:
        expand(JF_DIR + "depth_{meandepth}/{smp}.txt", smp = SMP_LIST, meandepth = ["{meandepth}"])
    output:
        MAT_DIR + "depth_{meandepth}/kmer-counts.tsv"
    log:
        MAT_DIR + "depth_{meandepth}/log-joinCounts.txt"
    shell:
        """
        echo -en "tag\\t" > {output}
        echo `ls {input} | sed 's/.txt//g' | sed 's|.*/sample|sample|g'` | sed 's/ /\\t/g' >> {output}
        apptainer exec -B "/data:/data" {KAMRAT_IMG} joinCounts -r 1 -a 1 {input} >> {output} 2> {log}
        """
