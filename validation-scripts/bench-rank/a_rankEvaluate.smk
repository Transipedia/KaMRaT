#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup snakemake --snakefile a_rankEvaluate.smk --cluster "qsub -q ssfa -l nodes=node29:ppn=1" --jobs 4 -p --latency-wait 60 --rerun-incomplete >> workflow_a_rankEvaluate.txt &

# Involved programs
COMPCODER = "/store/USERS/haoliang.xue/development/KaMRaT/paper-code/bench-rank/compcodeR.R"
RSCRIPT = "/home/haoliang.xue/.conda/envs/dekupl-hlx/bin/Rscript"
KAMRAT_IMG = "/home/haoliang.xue/tools/KaMRaT-2b3cf6e.sif"

# Outputs
RES_DIR = "/store/EQUIPES/SSFA/MEMBERS/haoliang.xue/KaMRaT-paper/benchmark-rank/"
IDX_DIR = "/data/work/I2BC/haoliang.xue/kamrat.idx/compcoder_simu/"
CPCDR_DIR = RES_DIR + "compcodeR_res/"
KAMRAT_DIR = RES_DIR + "kamrat_res/"

# Option lists
MODE_LIST = ["ft20000-es10-pout0", "ft20000-es1.5-pout0", "ft20000-es1.5-pout20", "ft200000-es1.5-pout20"]
MTHD_LIST = ["ttest.padj", "ttest.pi", "snr", "dids", "bayes", "lr"]

# ===== Workflow ===== #
rule all:
    input:
        expand(KAMRAT_DIR + "{mode}/ranked.{mthd}.tsv", mode = MODE_LIST, mthd = MTHD_LIST)
        
rule compcoder:
    output:
        CPCDR_DIR + "ft{nv}-es{es}-pout{pout}/countmat.tsv",
        CPCDR_DIR + "ft{nv}-es{es}-pout{pout}/smpinfo.tsv"
    log:
        CPCDR_DIR + "ft{nv}-es{es}-pout{pout}/log-compcoder.txt"
    params:
        nv = "{nv}",
        es = "{es}",
        pout = "{pout}",
        out = CPCDR_DIR + "ft{nv}-es{es}-pout{pout}/"
    shell:
        """
        {RSCRIPT} {COMPCODER} {params.nv} {params.es} {params.pout} {params.out} > {log}
        """

rule kamratIndex:
    input:
        CPCDR_DIR + "{mode}/countmat.tsv"
    output:
        IDX_DIR + "index-{mode}/idx-meta.bin",
        IDX_DIR + "index-{mode}/idx-pos.bin",
        IDX_DIR + "index-{mode}/idx-mat.bin"
    params:
        directory(IDX_DIR + "index-{mode}")
    log:
        IDX_DIR + "log-kamrat_index-{mode}.txt"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} -nfbase 30000000 2> {log}
        """

rule kamratRank:
    input:
        dsgn = CPCDR_DIR + "{mode}/smpinfo.tsv",
        idx1 = IDX_DIR + "index-{mode}/idx-meta.bin",
        idx2 = IDX_DIR + "index-{mode}/idx-pos.bin",
        idx3 = IDX_DIR + "index-{mode}/idx-mat.bin"
    output:
        KAMRAT_DIR + "{mode}/ranked.{mthd}.tsv"
    log:
        KAMRAT_DIR + "{mode}/log-ranked.{mthd}.txt"
    params:
        mthd = "{mthd}",
        idx = IDX_DIR + "index-{mode}"
    shell:
        """
        export SINGULARITY_BIND="/store:/store,/data:/data"
        singularity exec {KAMRAT_IMG} kamrat rank -idxdir {params.idx} -rankby {params.mthd} \
                                                  -design <(awk 'NR > 1 {{print $3"\\t"$1}}' {input.dsgn}) -outpath {output} -withcounts 2> {log}
        """
