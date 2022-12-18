# KaMRaT Validation Scripts

This folder collects the scripts that validate KaMRaT's functionality. These scripts reproduce the results in the demonstrating article (link to be added), and can be referred to as well for example usage.

The different goals for each script are listed below:

## Script Overview

```bash
bench-merge		                     # KaMRaT merge validation
|-- a_dataPrepare.smk                        #   Polyester simulation, k-mer count table preparation
|-- b_kamrat.smk                             #   KaMRaT index-merge, blastn on KaMRaT contigs
|-- c_rnaspades.smk                          #   rnaSPAdes assembling, blastn on rnaSPAdes contigs
|-- d1_misextension_ratio.R                  #   fig 2A
|-- d2_threshold-tuning.R                    #   fig 2B
|-- d3_cmp2spades.R                          #   fig 2C
bench-rank	                             # KaMRaT rank validation
|-- a_rankEvaluate.smk                       #   compcodeR simulation, KaMRaT index-rank
|-- b1_prcurves.R                            #   fig 3B
|-- b2_hclust.R                              #   fig 3C
|-- b3_upset.R                               #   fig S1
other-applications                           # KaMRaT application on LUADseo and PRADtcga datasets
|-- a_dataPrepare.smk                        #   adapter trimming, QC, and k-mer table preparation
|-- b_iMOKA.smk                              #   iMOKA application
|-- c_merge-rank.smk                         #   KaMRaT index-merge-rank (supervised)
|-- d_rank-merge.smk                         #   KaMRaT index-rank-merge (supervised)
|-- e_filter-merge.smk                       #   KaMRaT index-filter-merge
|-- e1_specctg-heatmap.R                     #   fig S5
|-- f_unsupv-rk-merge.smk                    #   KaMRaT index-rank-merge (unsupervised)
|-- g_corr-rk-merge.smk                      #   KaMRaT index-rank-merge (correlated features)
|-- g1_corr-heatmaps.R                       #   fig S4
|-- h1_barplots-efficiency-effectiveness.R   #   fig 4B, 4C, 4D
|-- h2_cmp_kamrat_imoka.R                    #   fig S2
|-- h3_cmp_kamrat-rank.R                     #   fig S3
```

## KaMRaT merge validation
### Required Data
[Gencode transcript reference (version 34)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz)

### Required R Library
[Polyester](https://bioconductor.org/packages/release/bioc/html/polyester.html)
