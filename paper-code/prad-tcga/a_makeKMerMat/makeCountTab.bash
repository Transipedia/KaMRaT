#!/bin/bash

/store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/makeKMerCountTab.sh --unlock --rerun-incomplete
nohup /store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/makeKMerCountTab.sh --configfile config.json --cluster "qsub -q ssfa -l nodes=1:ppn=3" --jobs 3 --rerun-incomplete --latency-wait 600 -p &
