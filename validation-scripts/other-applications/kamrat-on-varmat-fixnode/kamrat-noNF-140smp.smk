#!/bin/python3
# Author: Haoliang Xue

# nohup command:
#       nohup /home/haoliang.xue_ext/miniconda3/envs/kamrat-valid/bin/snakemake --snakefile kamrat-noNF-140smp.smk --cluster "qsub -q lowprio -l nodes=node40:ppn=1 -l mem=400G -l walltime=300:00:00 -m ea -M xhl-1993@hotmail.com" --jobs 6 -p --latency-wait 60 --rerun-incomplete >> workflow_kamrat-on-varmat2.txt &

# Involved programs
KAMRAT_IMG = "/home/haoliang.xue_ext/KaMRaT.sif"

# Parameters
SMPLIST_FILE = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-list.txt"
SMPCONDI_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/KaMRaT/validation-scripts/other-applications/sample-conditions/"
SMP_LIST = []
for f in ["LUADseo.tsv", "PRADtcga.tsv"]:
    with open(SMPCONDI_DIR + f) as scf:
        for row in scf:
            SMP_LIST.append(row.split("\t")[0])

# Inputs and outputs
RES_DIR = "/store/plateformes/CALCUL/SSFA_KaMRaT/Results/application-on-realdata/"
JF_DIR = RES_DIR + "jellyfish_res/"
MAT_DIR = RES_DIR + "joined_matrices/"
KAMRAT_DIR = "/data/work/I2BC/haoliang.xue/kamrat-new-res/Results/kamrat_on_varmat/"

rule all:
    input:
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/{f}.bin", \
               nsmp = [140], f = ["idx-meta", "idx-pos", "idx-mat"])

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

rule joinCounts_nsmp:
    input:
        meta = MAT_DIR + "all-smp.shuf.tsv",
        jf_res = expand(JF_DIR + "{smp}.txt", smp = SMP_LIST)
    output:
        meta = KAMRAT_DIR + "{n}-smp.tsv",
        tab = KAMRAT_DIR + "kmer-counts.{n}-smp.tsv.gz"
    log:
        KAMRAT_DIR + "log-joinCounts.{n}-smp.txt"
    params:
        n = "{n}",
        tab = KAMRAT_DIR + "kmer-counts.{n}-smp.tsv"
    threads: 1
    shell:
        """
        head -n{params.n} {input.meta} > {output.meta}
        awk 'BEGIN {{printf("tag")}} {{printf("\t%s", $1)}} END {{print ""}}' {output.meta} > {params.tab} # print the header row
        path_list=`awk -v jd="{JF_DIR}" '{{print jd$1".txt"}}' {output.meta}`
        apptainer exec -B "/store:/store" -B "/data:/data" {KAMRAT_IMG} joinCounts -r 1 -a 5 ${{path_list[@]}} >> {params.tab} 2> {log}
        gzip {params.tab}
        """

rule kamrat_index_nsmp:
    input:
        KAMRAT_DIR + "kmer-counts.{nsmp}-smp.tsv.gz"
    output:
        expand(KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/{f}.bin", \
               nsmp = ["{nsmp}"], f = ["idx-meta", "idx-pos", "idx-mat"])
    threads: 1
    log:
        KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/log-kamrat-index.txt"
    params:
        idx = KAMRAT_DIR + "apply-on-varmat/{nsmp}-smp/index/",
        nffile = KAMRAT_DIR + "140-NF.tsv"
    shell:
        """
        apptainer exec -B "/store:/store" -B "/data:/data" {KAMRAT_IMG} kamrat index -intab {input} -outdir {params.idx} \
                                                                                     -klen 31 -unstrand -nffile {params.nffile} &> {log}
        """

