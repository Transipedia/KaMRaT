#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile kamrat-index-merge.smk --cluster "qsub -q lowprio -l nodes=1:ppn=1 -l mem=200G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_kamrat-index-merge.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"
MKBLASTDB = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/makeblastdb"
BLASTN = "/home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/blastn"

# Parameters
SMPLIST_FILE = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-list.txt"
SMPCONDI_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-conditions/"
SMP_LIST = []
with open(SMPCONDI_DIR + "PRADtcga.tsv") as scf:
    for row in scf:
        SMP_LIST.append(row.split("\t")[0])

# Inputs and outputs
REF_FASTA = "/store/plateformes/CALCUL/SSFA_KaMRaT/Data/gencode.v34.transcripts.fa"
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
BLAST_DB = RES_DIR + "blast_db/" 
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "joined_matrices/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"

rule all:
    input:
        expand(KAMRAT_DIR + "merge-on-real/{dataset}/contig-alignment.{mode}.tsv",
               dataset = ["PRADtcga"],
               mode = ["mac", "pearson", "spearman", "none"])

rule makeblastdb:
    input:
        REF_FASTA
    output:
        BLAST_DB + "gencode.v34.transcripts.fa",
        BLAST_DB + "gencode.v34.transcripts.fa.nhr",
        BLAST_DB + "gencode.v34.transcripts.fa.nin",
        BLAST_DB + "gencode.v34.transcripts.fa.nog",
        BLAST_DB + "gencode.v34.transcripts.fa.nsd",
        BLAST_DB + "gencode.v34.transcripts.fa.nsi",
        BLAST_DB + "gencode.v34.transcripts.fa.nsq"
    threads: 1
    shell:
        """
        cp {REF_FASTA} {BLAST_DB}
        {MKBLASTDB} -in {input} -dbtype nucl -parse_seqids -out {BLAST_DB}gencode.v34.transcripts.fa
        """

rule joinCounts_datasets:
    input:
        smpcondi_file = SMPCONDI_DIR + "PRADtcga.tsv",
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        dsgn = MAT_DIR + "PRADtcga-smp.tsv",
        tab = MAT_DIR + "kmer-counts.PRADtcga-smp.tsv.gz"
    log:
        MAT_DIR + "log-joinCounts.PRADtcga-smp.txt"
    params:
        MAT_DIR + "kmer-counts.PRADtcga-smp.tsv"
    threads: 1
    shell:
        """
        cp {input.smpcondi_file} {output.dsgn}
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {output.dsgn} > {params} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {output.dsgn}`
        apptainer exec -B "/store:/store" {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {params} 2> {log}
        gzip {params}
        """

rule kamrat_index:
    input:
        MAT_DIR + "kmer-counts.{dataset}-smp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "merge-on-real/{dataset}/index/{f}.bin",\
               dataset = ["{dataset}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "merge-on-real/{dataset}/index/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "merge-on-real/{dataset}/index/"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                    -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_merge:
    input:
        files = expand(KAMRAT_DIR + "merge-on-real/{dataset}/index/{f}.bin",\
                       dataset = ["{dataset}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    output:
        KAMRAT_DIR + "merge-on-real/{dataset}/contig-counts.{mode}.tsv"
    threads: 1
    params:
        mode = "{mode}",
        idxdir = KAMRAT_DIR + "merge-on-real/{dataset}/index/"
    log:
        KAMRAT_DIR + "merge-on-real/{dataset}/log-kamrat-merge.txt"
    shell:
        """
        apptainer exec -B "/store:/store" {KAMRAT_IMG} kamrat merge -idxdir {params.idxdir} -overlap 30-15 \
                                                                    -interv {params.mode}:0.20 \
                                                                    -outpath {output} -withcounts rep &> {log}
        """

rule blastn:
    input:
        ctg_file = KAMRAT_DIR + "merge-on-real/{dataset}/contig-counts.{mode}.tsv",
        blast_db = BLAST_DB + "gencode.v34.transcripts.fa"
    output:
        KAMRAT_DIR + "merge-on-real/{dataset}/contig-alignment.{mode}.tsv"
    threads: 1
    params:
        KAMRAT_DIR + "merge-on-real/{dataset}/contig-seq.{mode}.fa"
    shell:
        """
        awk 'NR > 1 && $2 - 1 > 0 {{print ">ctg_"NR - 1"\\n"$1}}' {input.ctg_file} > {params}
        echo -e "qseqid\\tqlen\\tqstart\\tqend\\tslen\\tsstart\\tsend\\tlength\\tnident\\tpident\\tsseqid\\tqcovs\\tqacc" > {output}
        {BLASTN} -db {input.blast_db} -query {params} -max_hsps 1 -max_target_seqs 1 -evalue 0.001 \
                 -outfmt "6 qseqid qlen qstart qend slen sstart send length nident pident sseqid qcovs qacc" -num_threads 1 -dust no >> {output}
        """
