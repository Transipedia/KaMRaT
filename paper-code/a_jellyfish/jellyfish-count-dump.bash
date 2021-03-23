#!/bin/bash

/store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/jellyfish-count-dump.sh --unlock --rerun-incomplete
nohup /store/USERS/haoliang.xue/development/KaMRaT/related-tools/prepare_kmer_table/jellyfish-count-dump.sh --configfile config.json --cluster "qsub -q common -l nodes=1:ppn=3" --jobs 3 --rerun-incomplete --latency-wait 600 -p &
