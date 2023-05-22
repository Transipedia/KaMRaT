#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile e_filter-merge.smk --cluster "qsub -q common -l nodes=node27:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr" --jobs 1 -p --latency-wait 60 --rerun-incomplete >> workflow_e_filter-merge.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"
BLASTN = "/usr/bin/blastn"

configfile: "luad-config.json"

# Outputs
RES_DIR = config["res_dir"]
MAT_DIR = RES_DIR + "matrices/"
DSGN_PATH = MAT_DIR + "all-samples.tsv"
FILTER_MERGE_DIR = RES_DIR + "kamrat_res/all-samples/filter-merge/"
IDX_DIR = "/data/work/I2BC/haoliang.xue/kamrat.idx/LUADseo/all-samples/index-raw/"
BLAST_DB = "/store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa"

rule all:
    input:
        expand(FILTER_MERGE_DIR + "spec-contig-counts.tsv")

rule index:
    input:
        MAT_DIR + "kmer-counts.allsmp.tsv"
    output:
        IDX_DIR + "idx-mat.bin",
        IDX_DIR + "idx-meta.bin",
        IDX_DIR + "idx-pos.bin"
    threads: 1
    log:
        IDX_DIR + "log-kamratIndex.allsmp.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat index -intab {input} -outdir {IDX_DIR} -klen 31 -unstrand &> {log} # no normalization
        """

rule filter_merge:
    input:
        IDX_DIR + "idx-mat.bin",
        IDX_DIR + "idx-meta.bin",
        IDX_DIR + "idx-pos.bin"
    output:
        ft = FILTER_MERGE_DIR + "spec-kmers.bin",
        mg = FILTER_MERGE_DIR + "spec-contig-counts.tsv"
    threads: 1
    log:
        FILTER_MERGE_DIR + "log-filter-merge.allsmp.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat filter -idxdir {IDX_DIR} -design <(sed 's/normal/DOWN/g' {DSGN_PATH} | sed 's/tumor/UP/g') -upmin 10:12 -downmax 0:77 -outpath {output.ft} &> {log}
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {IDX_DIR} -overlap 30-15 -with {output.ft} -interv pearson:0.20 -outpath {output.mg} -withcounts mean &>> {log}
        """

