# Companion scripts to KaMRaT

Here we provide several companion scripts for a more friendly usage of KaMRaT.

## dekupl-joinCounts

This is a submodule cloned from [here](https://github.com/Transipedia/dekupl-joinCounts/) which was initially developped for the [DE-kupl software](https://doi.org/10.1186/s13059-017-1372-2).

This submodule merges Jellyfish count results of multiple samples into a single k-mer count matrix, and is called by the Snakefile in make-matrix folder (see below). Please refer to dekupl-joinCounts GitHub page for detailed information and cite [DE-kupl article](https://doi.org/10.1186/s13059-017-1372-2) for the usage of it.

## make-matrix

This folder includes a Snakefile summarising workflow to produce k-mer count matrix from fastq files.

A perl script `revCompFastq.pl` from [DE-kupl GitHub](https://github.com/Transipedia/dekupl-run) is also included, to generate reverse complement sequence for reverse reads.

To run the snakemake workflow, all the fastq files are required to be wrapped in a same folder, with suffixes in the same pattern. For example, in this show case, all the fastq files for the two toy samples are wrapped in the folder `toyroom/data/fastq_dir/`, with the shared suffixes pattern as `.R1.fastq.gz` and `.R2.fastq.gz`.

Users are expected to prepare and provide a `config.json` file for the `snakemake` workflow. An example can be found in this folder named `toy-config.json`, with the keys being explained as below:

- samples_tsv: a file indicating samples to be analysed. Multi-column tab-separated table is allowed, but a header line with the first field being "sample" is mandatory.
- lib_type: sequencing strandedness, can be one of "rf", "fr" or unstranded.
- fastq_dir: the folder that wraps fastq files related to the samples provided in `samples_tsv`.
- r1_suffix: suffix patterns of the first read files.
- r2_suffix: suffix patterns of the second read files.
- output_dir: output folder.
- kmer_length: k-mer length.
- min_rec: requirement of minimum recurrent sample number to report the k-mer.
- min_rec_abd: requirement of minimum abundance threshold to report k-mer occurrence in one sample.
- n_cores: number of CPU cores to be used for computation.

To launch the workflow, please run:

```bash
apptainer exec -B /in_dir/:/sif_data/ -B /out_dir/:/sif_out/ -B $PWD:/sif_pwd/ KaMRaT.sif \
               snakemake -s /usr/KaMRaT/related-tools/make-matrix/Snakefile \
                         --configfile /sif_pwd/toy-config.json --cores 1
```

A toy use case with examples of `config.json` and `sample_tsv` files can be found in the `toyroom/usecases/MakeTab.usecases/` folder of this GitHub repository.

## splitCV.bash

A bash script shuffling and splitting feature count matrix for n-fold CV analysis.

## selectSubmat.bash

A bash script for retrieving specific samples (in columns) from a given matrix.

## zero-counter.bash

An one-line command in AWK for counting number of zeros in the given matrix.

## downstream_analysis

Scripts for downstream analyses, e.g., classifier construction.
