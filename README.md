# KaMRaT

-----
KaMRaT is a C++ tool for finding substrings with interesting properties in large NGS datasets. 

KaMRaT requires a k-mer count matrix extracted from the NGS files (e.g. with Jellyfish), and potential design table or FASTA file (see below for more information). 

KaMRaT then provides a set of tools for reducing the k-mer matrix and extending k-mers to longer contigs. The main modules are:

- kamrat index: index feature* count table on disk
- kamrat filter: remove/retain features* by expression level 
- kamrat mask: select and suppress k-mers matching given fasta sequences
- kamrat merge: merge k-mers into contigs
- kamrat score: score features* by classification performance, statistical significance, correlation, or variability 
- kamrat query: estimate count vectors of given list of contigs

  Note: \*	features can be not only k-mers or k-mer contigs, but also general features such as genes or transcripts.

KaMRaT means "k-mer Matrix Reduction Toolkit", or "k-mer Matrix, Really Tremendous !".

KaMRaT has been published on [*Bioinformatics*](https://academic.oup.com/bioinformatics/article/40/3/btae090/7623008).

**Before using KaMRaT, please refer to our [License](License.md).**

-----

**Note: this page contains a TL;DR version of the software installation and usages. Please refer to our [Wiki page](https://github.com/Transipedia/KaMRaT/wiki) for detailed information.**

## Software structure
KaMRaT is designed as a flexible toolkit for combinations of different operations. The workflow always starts from index, then other operations can be combined according to user's analysis design, including:

- a single operation of filter, mask, merge, score or query;
- filter-merge or filter-score;
- mask-merge or mask-score;
- merge-score;
- score-merge.

![workflow](./docs/workflow.png)

### Tips for the choice of k-mer length
Currently, KaMRaT only accepts k-mers no longer than 32nt, since the k-mers are coded in an uint64 variable.

Besides, we recommend users to choose k as an odd number, to avoid confounding one k-mer with its reverse complement counterpart in unstranded data. For example, in the situation k=6, 6-mers such as `AAATTT` lose their information of strandedness.


## Current Release KaMRaT v1.2
The current release of KaMRaT is v1.2. Compared to its previous release, v1.1, it introduces several new characteristics:
- Each module of `filter`, `merge` and `score` now supports outputting a fasta file containing selected/merged sequences.
- The output tables of all modules `filter`, `mask`, `merge`, `score` and `query` now output the count table with values being rounted to the nearest integers. The decimal values can be output by setting `-counts` argument.
- The `index` module now allows a normalisation factor base too large or too small, instead of throwing an exception, it puts a warning.
- The `mask` module now allows simultaneously indicating sequences to select and suppress.
- Updated license.
- Bugfix: the `merge` module was mistakenly output '<' as fasta header leading character, corrected to '>'.

For full release notes, please refer to our [wiki page](https://github.com/Transipedia/KaMRaT/wiki/General-Descriptions).


## Installation
It's highly recommended to directly use KaMRaT within `apptainer`/`singularity` container for users at any level unless the task involves in software development, because:

1. Building the container is simple and only requires the dependency of `apptainer`/`singularity` which should be pre-installed on most of HPC clusters.
2. Usage of KaMRaT in container ensures better reproducibility of results.
3. Several companion scripts can be easily run within the container, e.g., for input k-mer matrix construction.

Pre-built image is provided on `DockerHub`. To build in local, please run:

```bash
apptainer build KaMRaT.sif docker://xuehl/kamrat:latest
    # alternatively to build singularity image: simply replace "apptainer" to "singularity"
```

Please refer to our [Wiki page](https://github.com/Transipedia/KaMRaT/wiki/Software-Installation) for:
- more detailed information on KaMRaT usage within `apptainer`/`singularity` container;
- alternatively way of installation by building the software from source.

## Usage
KaMRaT can be run generally in the fashion:

``` bash
apptainer exec -B /bind_src:/bind_des kamrat <CMD> [options] /path/from/{bind_des}/to/input/kmer/table 
    # <CMD> can be one of index, filter, mask, merge, score, query
    # replace "apptainer" to "singularity" when KaMRaT is built by singularity
```

The top-level helper is reachable by:

``` bash
apptainer exec kamrat
```

Helpers of specific KaMRaT modules are accessible via:

``` bash
apptainer exec kamrat <CMD>
    # <CMD> can be one of index, filter, mask, merge, score, query
```

Please refer to our [Wiki page](https://github.com/Transipedia/KaMRaT/wiki/Software-Usage) for detailed software usage description.

Demostrations of example usecases can be found [here](https://github.com/Transipedia/KaMRaT/wiki/Workflow-Demos).

## Software/Library Citations

Armadillo:

+ Conrad Sanderson and Ryan Curtin. Armadillo: a template-based C++ library for linear algebra. Journal of Open Source Software, Vol. 1, pp. 26, 2016.
+ Conrad Sanderson and Ryan Curtin. A User-Friendly Hybrid Sparse Matrix Class in C++. Lecture Notes in Computer Science (LNCS), Vol. 10931, pp. 422-430, 2018.

[Boost C++ Library](https://www.boost.org/)

DE-kupl: Audoux, J., Philippe, N., Chikhi, R. et al. DE-kupl: exhaustive capture of biological variation in RNA-seq data through k-mer decomposition. Genome Biol 18, 243, 2017.

MLPack: R.R. Curtin, M. Edel, M. Lozhnikov, Y. Mentekidis, S. Ghaisas, S. Zhang. mlpack 3: a fast, flexible machine learning library. Journal of Open Source Software 3:26, 2018.
