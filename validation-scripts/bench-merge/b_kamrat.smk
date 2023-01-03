#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile b_kamrat.smk --cluster "qsub -q common -l nodes=node15 -l mem=50G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 2 -p --latency-wait 60 --rerun-incomplete >> workflow_b_kamrat.txt &

# Involved programs
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"

# Inputs
REF_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Data/bench-merge/"

# Outputs
RES_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge/"
MAT_DIR = RES_DIR + "matrices/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"

MEAN_DEPTH_LIST = [0.1, 0.2, 0.5, 1, 2, 5, 10]
INTERV_LIST = ["pearson", "spearman", "mac"]
THRES_LIST = [str(i / 10) for i in range(1, 11)]

# ===== Workflow ===== #
rule all:
    input:
        REF_DIR + "gencode.v34.transcripts.fa.nhr",
        expand(KAMRAT_DIR + "depth_{mean_depth}/ctg-aligned.none.tsv", 
               mean_depth=MEAN_DEPTH_LIST),
        expand(KAMRAT_DIR + "depth_{mean_depth}/ctg-aligned.{interv}_{thres}.tsv",
               mean_depth=MEAN_DEPTH_LIST, interv=INTERV_LIST, thres=THRES_LIST)

rule makeblastdb:
    input:
        REF_DIR + "gencode.v34.transcripts.fa"
    output:
        REF_DIR + "gencode.v34.transcripts.fa.nhr",
        REF_DIR + "gencode.v34.transcripts.fa.nin",
        REF_DIR + "gencode.v34.transcripts.fa.nog",
        REF_DIR + "gencode.v34.transcripts.fa.nsd",
        REF_DIR + "gencode.v34.transcripts.fa.nsi",
        REF_DIR + "gencode.v34.transcripts.fa.nsq"
    threads: 1
    shell:
        """
        {MKBLASTDB} -in {input} -dbtype nucl -parse_seqids
        """

rule kamrat_index:
    input:
        MAT_DIR + "depth_{mean_depth}/kmer-counts.tsv"
    output:
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-meta.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-pos.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-mat.bin"
    params:
        KAMRAT_DIR + "depth_{mean_depth}/index/"
    log:
        KAMRAT_DIR + "depth_{mean_depth}/index/log-kamrat_index.txt"
    threads: 1
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat index -intab {input} \
                                                                  -outdir {params} \
                                                                  -klen 31 -unstrand -nfbase 1000000000 &> {log}
        """

rule kamrat_merge_none:
    input:
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-meta.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-pos.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-mat.bin"
    output:
        KAMRAT_DIR + "depth_{mean_depth}/merged-counts.none.tsv"
    params:
        KAMRAT_DIR + "depth_{mean_depth}/index/"
    threads: 1
    log:
        KAMRAT_DIR + "depth_{mean_depth}/log-kamrat-merge-none.txt"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params} \
                                                                  -overlap 30-15 -interv none \
                                                                  -outpath {output} -withcounts rep &> {log}
        """

rule kamrat_merge_interv:
    input:
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-meta.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-pos.bin",
        KAMRAT_DIR + "depth_{mean_depth}/index/idx-mat.bin"
    output:
        KAMRAT_DIR + "depth_{mean_depth}/merged-counts.{interv}_{thres}.tsv"
    log:
        KAMRAT_DIR + "depth_{mean_depth}/log-kamrat-merge-{interv}_{thres}.txt"
    threads: 1
    params:
        idxdir = KAMRAT_DIR + "depth_{mean_depth}/index/",
        interv = "{interv}",
        thres = "{thres}"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params.idxdir} \
                                                                  -overlap 30-15 -interv {params.interv}:{params.thres} \
                                                                  -outpath {output} -withcounts rep &> {log}
        """

rule blastn:
    input:
        ctg_file = KAMRAT_DIR + "depth_{mean_depth}/merged-counts.{mode}.tsv",
        blast_db = REF_DIR + "gencode.v34.transcripts.fa"
    output:
        KAMRAT_DIR + "depth_{mean_depth}/ctg-aligned.{mode}.tsv"
    threads: 1
    params:
        KAMRAT_DIR + "depth_{mean_depth}/ctg-seq.{mode}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input.ctg_file} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid" -num_threads 1 -dust no >> {output}
        """
