# Companion scripts to KaMRaT

Here we provide several companion scripts for a more friendly usage of KaMRaT.

## dekupl-joinCounts

This is a submodule cloned from [here](https://github.com/Transipedia/dekupl-joinCounts/) which was initially developped for the [DE-kupl software](https://doi.org/10.1186/s13059-017-1372-2).

This submodule merges Jellyfish count results of multiple samples into a single k-mer count matrix, and is called by the Snakefile in make-matrix folder (see below). Please refer to dekupl-joinCounts GitHub page for detailed information and cite [DE-kupl article](https://doi.org/10.1186/s13059-017-1372-2) for the usage of it.

## make-matrix

This folder includes a Snakefile summarising workflow to produce k-mer count matrix from fastq files.

A perl script `revCompFastq.pl` from [DE-kupl GitHub](https://github.com/Transipedia/dekupl-run) is also included, to generate reverse complement sequence for reverse reads.

To run the snakemake workflow, please refer to the `toyroom/usecases/MakeTab.usecases/` folder of this GitHub repository for a toy use case with examples of `config.json` and `sample_tsv` files.

## splitCV.bash

A bash script shuffling and splitting feature count matrix for n-fold CV analysis.

## selectSubmat.bash

A bash script for retrieving specific samples (in columns) from a given matrix.

## zero-counter.bash

An one-line command in AWK for counting number of zeros in the given matrix.

## downstream_analysis

Scripts for downstream analyses, e.g., classifier construction.
