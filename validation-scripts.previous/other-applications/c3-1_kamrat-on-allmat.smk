#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile c3-1_kamrat-on-allmat.smk --cluster "qsub -q lowprio -l nodes=1:ppn=1 -l mem=350G -l walltime=300:00:00 -m ea -M xue.haoliang@outlook.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_c3-1_kamrat-on-allmat.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Inputs and outputs
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
MAT_DIR = RES_DIR + "joined_matrices/"
KAMRAT_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/kamrat_on_varmat/"
NSMP = ["200"]

rule all:
    input:
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/kamrat-index-merge.pearson.tsv", nsmp = NSMP),
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/kamrat-index-rank.snr.tsv", nsmp = NSMP)

rule shuffle_samples:
    input:
        MAT_DIR + "all-smp.tsv"
    output:
        MAT_DIR + "all-smp.shuf.tsv"
    threads: 1
    log:
        MAT_DIR + "log-shuffle-samples.txt"
    shell:
        """
        get_seeded_random() {{
            seed="$1"
            openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
        }}
        sed 's/normal/A/g' {input} | sed 's/tumor/B/g' | sed 's/No/A/g' | sed 's/Yes/B/g' | shuf --random-source=<(get_seeded_random 91400) > {output} 2> {log}
        """

rule subset_mat:
    input:
        meta = MAT_DIR + "all-smp.shuf.tsv",
        tab = MAT_DIR + "kmer-counts.all-smp.tsv.gz"
    output:
        meta = KAMRAT_DIR + "{nsmp}-smp.tsv",
        tab = KAMRAT_DIR + "kmer-counts.{nsmp}-smp.tsv.gz"
    params:
        nsmp = "{nsmp}",
        tab = KAMRAT_DIR + "kmer-counts.{nsmp}-smp.tsv"
    threads: 1
    log:
        KAMRAT_DIR + "log-subset-mat.{nsmp}-smp.txt"
    shell:
        """
        head -n{params.nsmp} {input.meta} > {output.meta}
        col_list_str=`awk '{{print $1}}' {output.meta}`
        zcat {input.tab} | awk -v col2print="$col_list_str" 'BEGIN {{nc2p = split(col2print, c2p)}} NR == 1 {{for (i = 1; i <= NF; ++i) hname[$i]=i}} {{printf("%s", $1); for (j = 1; j <= nc2p; ++j) if (c2p[j] in hname) printf("\\t%s", $(hname[c2p[j]])); print ""}}' > {params.tab}
        gzip {params.tab}
        """

rule kamrat_index:
    input:
        KAMRAT_DIR + "kmer-counts.{nsmp}-smp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/{f}.bin", \
               nsmp = ["{nsmp}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/log-kamrat-index.txt"
    params:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/"
    shell:
        """
        apptainer exec -B "/store:/store" -B "/data:/data" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params} \
                                                                                     -klen 31 -unstrand -nfbase 2000000000 &> {log}
        """

rule kamrat_merge:
    input:
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/{f}.bin",
               nsmp = ["{nsmp}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    output:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/kamrat-index-merge.{mode}.tsv"
    threads: 1
    resources:
        mem_mb = 350 * 1024 # 350G memory
    params:
        idxdir = KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/",
        mode = "{mode}"
    log:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/log-kamrat-merge-{mode}.txt"
    shell:
        """
        apptainer exec -B "/store:/store" -B "/data:/data" {KAMRAT_IMG} kamrat merge -idxdir {params.idxdir} -overlap 30-15 \
                                                                                     -interv {params.mode}:0.20 -outpath {output} \
                                                                                     -withcounts mean &> {log}
        """

rule kamrat_rank:
    input:
        meta = KAMRAT_DIR + "{nsmp}-smp.tsv",
        files = expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/{f}.bin",
                       nsmp = ["{nsmp}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    output:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/kamrat-index-rank.{mode}.tsv"
    threads: 1
    resources:
        mem_mb = 350 * 1024 # 350G memory
    params:
        idxdir = KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/",
        mode = "{mode}"
    log:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/log-kamrat-rank-{mode}.txt"
    shell:
        """
        apptainer exec -B "/store:/store" -B "/data:/data" {KAMRAT_IMG} kamrat rank -idxdir {params.idxdir} -rankby {params.mode} \
                                                                                    -design {input.meta} -outpath {output} \
                                                                                    -withcounts &> {log}
        """
