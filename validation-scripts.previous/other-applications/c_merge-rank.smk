#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile c_merge-rank.smk --cluster "qsub -l nodes=node29:ppn=1 -m ea -M haoliang.xue@i2bc.paris-saclay.fr" --jobs 3 -p --latency-wait 60 --rerun-incomplete >> workflow_c_merge-rank_luad.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"
IMOKA_IMG = "/home/haoliang.xue/tools/iMOKA-1.1.img"
RSCRIPT = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript"
TRTAB4IMOKA = "/store/USERS/haoliang.xue/development/KaMRaT/paper-code/real-data/reformat4iMOKA.R"

configfile: "luad-config.json"
# configfile: "prad-config.json"

# Outputs
RES_DIR = config["res_dir"]
MAT_DIR = RES_DIR + "matrices/"
TMP_DIR = config["tmp_dir"]
KAMRAT_DIR = RES_DIR + "kamrat_res/"
TOP_NUM = config["top_num"]

MTHD_LIST = ["ttest.padj", "ttest.pi", "snr", "lr", "dids", "bayes"]
OTHER_MODE_LIST = ["none", "mac:0.3", "spearman:0.2"]

rule all:
    input:
        expand(KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/prediction_test.txt", num = range(5), mthd = MTHD_LIST),
        expand(TMP_DIR + "CV{num}/merged-contigs-{mode}.bin", num = range(5), mode = OTHER_MODE_LIST)

rule kamrat_index:
    input:
        smpinfo = MAT_DIR + "sampleshuf.{set}{num}.tsv",
        kmercounts = MAT_DIR + "kmer-counts.{set}{num}.tsv"
    output:
        directory(TMP_DIR + "CV{num}/index-{set}")
    threads: 1
    log:
        KAMRAT_DIR + "CV{num}/merge-rank/log-kamrat-index-{set}{num}.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        mkdir -p {output}
        singularity exec {KAMRAT_IMG} kamrat index -intab {input.kmercounts} -outdir {output} -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_merge:
    input:
        TMP_DIR + "CV{num}/index-train"
    output:
        mg = TMP_DIR + "CV{num}/merged-contigs.bin"
    threads: 1
    resources:
        mem_mb = 120 * 1024 # 120G memory
    log:
        KAMRAT_DIR + "CV{num}/merge-rank/log-kamrat-merge.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {input} -overlap 30-15 -interv pearson:0.20 -outpath {output.mg} &> {log}
        """

rule kamrat_merge_others:
    input:
        TMP_DIR + "CV{num}/index-train"
    output:
        mg = TMP_DIR + "CV{num}/merged-contigs-{mode}.bin"
    threads: 1
    resources:
        mem_mb = 120 * 1024 # 120G memory
    log:
        KAMRAT_DIR + "CV{num}/merge-rank/log-kamrat-merge-{mode}.txt"
    params:
        "{mode}"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat merge -idxdir {input} -overlap 30-15 -interv {params} -outpath {output.mg} &> {log}
        """

rule kamrat_rank:
    input:
        smpinfo = MAT_DIR + "sampleshuf.train{num}.tsv",
        idx = TMP_DIR + "CV{num}/index-train",
        mg = TMP_DIR + "CV{num}/merged-contigs.bin"
    output:
        rk = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/top-merged-contig-counts.tsv"
    params:
        mthd = "{mthd}"
    threads: 1
    resources:
        mem_mb = 120 * 1024 # 120G memory
    log:
        KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/log-kamrat-merge-rank.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat rank -idxdir {input.idx} -rankby {params} -design {input.smpinfo} -with {input.mg}:mean \
                                                  -seltop {TOP_NUM} -outpath {output.rk} -withcounts &> {log}
        """

rule kamrat_build_and_test:
    input:
        rk = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/top-merged-contig-counts.tsv",
        idx_test = TMP_DIR + "CV{num}/index-test",
        smpinfo_train = MAT_DIR + "sampleshuf.train{num}.tsv",
        smpinfo_test = MAT_DIR + "sampleshuf.test{num}.tsv"
    output:
        tab4imoka_train = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/top-merged-contig-counts.4imoka.tsv",
        rf = directory(KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/randomforest_models/"),
        qres = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/query-res-intest.tsv",
        tab4imoka_test = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/query-res-intest.4imoka.tsv",
        pred = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/prediction_test.txt"
    params:
        mthd = "{mthd}",
        rf_out = KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/randomforest"
    threads: 1
    log:
        KAMRAT_DIR + "CV{num}/merge-rank/{mthd}/log-kamrat-build-and-test.txt"
    shell:
        """
        export IMOKA_MAX_MEM_GB=100
        export SINGULARITY_BINDPATH="/store:/store,/data:/data"
        export OMP_NUM_THREADS=4
        {RSCRIPT} {TRTAB4IMOKA} {input.rk} {input.smpinfo_train} {output.tab4imoka_train} &> {log}
        singularity exec {IMOKA_IMG} random_forest.py -m 100 {output.tab4imoka_train} {params.rf_out} &>> {log}
        awk '{{print ">contig"NR"\\n"$1}}' {output.rf}0_RF.pickle.features > {output.rf}0_RF.pickle.features.fa
        singularity exec {KAMRAT_IMG} kamrat query -idxdir {input.idx_test} -fasta {output.rf}0_RF.pickle.features.fa -toquery median -withabsent -outpath {output.qres} &>> {log}
        {RSCRIPT} {TRTAB4IMOKA} {output.qres} {input.smpinfo_test} {output.tab4imoka_test} &>> {log}
        singularity exec {IMOKA_IMG} predict.py {output.tab4imoka_test} {output.rf}0_RF.pickle {output.pred} &>> {log}
        """
