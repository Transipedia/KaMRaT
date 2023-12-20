# Making k-mer count matrix from fastq files

The singularity image of KaMRaT integrates a `snakemake` workflow to produce k-mer count matrix from fastq files. Here a toy example is provided as an instruction.

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

An example value for the `samples_tsv` key is provided in `sample.tsv` file.

To launch the workflow, please run:

```bash
singularity exec -B /in_dir/:/sif_data/ -B /out_dir/:/sif_out/ -B $PWD:/sif_pwd/ KaMRaT.sif \
                 snakemake -s /usr/KaMRaT/related-tools/make-matrix/Snakefile \
                           --configfile /sif_pwd/toy-config.json --cores 1
```
