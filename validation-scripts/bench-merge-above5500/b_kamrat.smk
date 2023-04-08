#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile b_kamrat.smk --cluster "qsub -q common -l nodes=1:ppn=3 -l mem=100G -l walltime=300:00:00" --jobs 3 -p --latency-wait 60 --rerun-incomplete >> workflow_b_kamrat.txt &

# Involved programs
RSCRIPT = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"

# Inputs and outputs
RES_DIR = f"/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/bench-merge-above5500/"
SMPREF_DIR = RES_DIR + "smpref_res/"
BLAST_DB = RES_DIR + "blast_db/"
MAT_DIR = RES_DIR + "matrices/"

# Parameters
MEAN_DEPTH_LIST = [0.05, 0.2, 0.4, 0.6, 0.8, 1]
INTERV_LIST = ["pearson", "spearman", "mac"]
THRES_LIST = [str(i / 10) for i in range(1, 11, 2)]


# ===== Workflow ===== #
rule all:
    input:
        expand(RES_DIR + "kamrat_res_err-free_1-1/depth_{mean_depth}/ctg-aligned.none.tsv", 
               mean_depth=MEAN_DEPTH_LIST),
        expand(RES_DIR + "kamrat_res_err-free_1-1/depth_{mean_depth}/ctg-aligned.{interv}_0.2.tsv",
               mean_depth=MEAN_DEPTH_LIST, interv=INTERV_LIST),
        expand(RES_DIR + "kamrat_res_err-free_1-1/depth_0.05/ctg-aligned.{interv}_{thres}.tsv",
               interv=INTERV_LIST, thres=THRES_LIST)

rule makeblastdb:
    output:
        BLAST_DB + "merged-ref.fa",
        BLAST_DB + "merged-ref.fa.nhr",
        BLAST_DB + "merged-ref.fa.nin",
        BLAST_DB + "merged-ref.fa.nog",
        BLAST_DB + "merged-ref.fa.nsd",
        BLAST_DB + "merged-ref.fa.nsi",
        BLAST_DB + "merged-ref.fa.nsq"
    threads: 1
    shell:
        """
        for f in {SMPREF_DIR}*.fa
        do
            f1=$(basename $f ".fa")
            sed "s/>ENST/>$f1|ENST/g" $f >> {BLAST_DB}merged-ref.fa
        done
        {MKBLASTDB} -in {BLAST_DB}merged-ref.fa -dbtype nucl -parse_seqids -out {BLAST_DB}merged-ref.fa
        """

rule kamrat_index:
    input:
        MAT_DIR + "{errmode}/depth_{mean_depth}/kmer-counts-{min_rec}-{min_abd}.tsv"
    output:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-meta.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-pos.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-mat.bin"
    params:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/"
    log:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/log-kamrat_index.txt"
    threads: 1
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat index -intab {input} \
                                                                  -outdir {params} \
                                                                  -klen 31 -unstrand -nfbase 500000000 &> {log}
        """

rule kamrat_merge_none:
    input:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-meta.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-pos.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-mat.bin"
    output:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/merged-counts.none.tsv"
    params:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/"
    threads: 1
    log:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/log-kamrat-merge-none.txt"
    shell:
        """
        apptainer exec -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params} \
                                                                  -overlap 30-15 -interv none \
                                                                  -outpath {output} -withcounts rep &> {log}
        """

rule kamrat_merge_interv:
    input:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-meta.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-pos.bin",
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/idx-mat.bin"
    output:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/merged-counts.{interv}_{thres}.tsv"
    log:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/log-kamrat-merge-{interv}_{thres}.txt"
    threads: 1
    params:
        idxdir = RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/index/",
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
        ctg_file = RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/merged-counts.{mode}.tsv",
        blast_db = BLAST_DB + "merged-ref.fa"
    output:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/ctg-aligned.{mode}.tsv"
    threads: 3
    params:
        RES_DIR + "kamrat_res_{errmode}_{min_rec}-{min_abd}/depth_{mean_depth}/ctg-seq.{mode}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input.ctg_file} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tqcovs\\tsseqid" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -max_hsps 1 -max_target_seqs 1 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident qcovs sseqid" -num_threads 3 -dust no >> {output}
        """
