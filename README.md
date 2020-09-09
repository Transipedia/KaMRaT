# KaMRaT

KaMRaT means "k-mer Matrix Reduction Toolbox", or "k-mer Matrix, Really Tremendous !".

The toolbox contains for now modules below:

- kamratMerge which takes a k-mer count table as input, and merge k-mers into longer contigs according to their overlap
- kamratReduce which takes a k-mer count table as input, and evaluates the performance of each k-mer according to different metrics

## Prerequest

The kamratReduce module is dependent on MLPack library. It could be installed from [conda cloud](https://anaconda.org/conda-forge/mlpack).

After installation, please add the following line into your ```.bashrc``` file in the ```home/``` directory:

```bash
export LD_LIBRARY_PATH="/path_to_conda_env/mlpack/lib:$LD_LIBRARY_PATH"
```

For compiling the source, you can use ```compile.bash``` in root directory of KaMRaT:

```bash
bash compile.bash
```

And the executable files are in ```bin/``` directory.

## Convention

### File of Sample Info

The sample-info file is indicated by the option ```-d```. This file aims to indicate which columns in the k-mer count matrix should be considered as sample columns. Please do not add any header line in the file, since the column contents are already defined as below.

If provided, this file may contains one or two columns:

- if the file contains only one column, it indicates sample names, and all samples are considered as in same condition
- if the file contains two columns, the first column correspond to sample names, and the second conrresponds to conditions

If not provided, all columns apart from the first one in the k-mer count matrix are considered as samples.
