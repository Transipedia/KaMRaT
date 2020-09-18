# KaMRaT

k-mers are the substrings of a fixed length $k$ along biological sequences.

Transcriptomics analysis based on k-mer count possesses good potential for capturing signals at nucleotide solution. But a current challenge is that the number of k-mer features are often too many to have an effective analysis.

KaMRaT provides a set of tools for k-mer analysis, for reducing number of k-mer features in consideration. The name KaMRaT means "k-mer Matrix Reduction Toolkit", or "k-mer Matrix, Really Tremendous !".

## Before KaMRaT Installation

### Dependencies

KaMRaT is dependent on the following software/libraires:

- [Git](https://git-scm.com/)
- If running KaMRaT in Singularity
  - [Singularity](https://sylabs.io/docs/)
- If running KaMRaT without container
  - [MLPack 3.3.2](https://github.com/mlpack/mlpack/releases/tag/3.3.2)
  - [Boost-iostreams](https://www.boost.org/doc/libs/1_74_0/libs/iostreams/doc/index.html)

The simplest way for KaMRaT installation is by using Singularity container which can be installed with [the guidance here](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps). In this case, you can jump to [KaMRaT Installation](##KaMRaT-Installation)

If you prefer to run KaMRaT without any container, you can follow the guidance below.

### MLPack Installation

MLPack can be installed on [Linux/Mac](https://mlpack.org/doc/mlpack-3.3.2/doxygen/build.html), [Windows](https://mlpack.org/doc/mlpack-3.3.2/doxygen/build_windows.html), or via [conda](https://anaconda.org/conda-forge/mlpack) by following the corresponding links.

If you are installing MLPack with conda, please add the following line into your ```.bashrc``` file in the ```home/``` directory before compiling KaMRaT:

```bash
export LD_LIBRARY_PATH="/path_to_conda_env/mlpack/lib:$LD_LIBRARY_PATH"
```

### Boost-iostreams Library Installation

Boost-iostreams library can be installed by the following command:

```bash
apt-get install libboost-iostreams-dev
```

## KaMRaT Installation

### Installation Option 1: inside singularity

```bash
git clone https://github.com/Transipedia/KaMRaT.git
cd KaMRaT
mkdir bin
singularity build --fakeroot bin/kamrat.simg kamrat.simg.def
```

A sigularity image is created in the ```bin/``` directory.

### Installation Option 2: compiling from source

Compiling command:

```bash
git clone https://github.com/Transipedia/KaMRaT.git
cd KaMRaT
```

Then, if you installed MLPack library with conda:

```bash
bash compile.bash /path_to_MLPack_conda_environment
```

Otherwise, if you installed MLPack without conda:

```bash
bash compile.bash
```

An executable binary file is available in the ```bin/``` directory.

## Conventions Before KaMRaT Execution

### File of Sample Info

The sample-info file is indicated by the option ```-smp-info```. This file aims to indicate which columns in the k-mer count matrix should be considered as sample columns. Please do not put any header line in the file, since the columns are already defined as below.

If provided, this file may contains one or two columns:

- If the file contains only one column, it indicates sample names, and all samples are considered as the same condition
- If the file contains two columns, the first column corresponds to sample names, and the second conrresponds to conditions

If not provided, all columns apart from the first one in the matrix are considered as samples.

### Metrics for naïve Bayes Evaluation in KaMRaT rank

KaMRaT supports ranking and reducing features by using naïve Bayes model with F1-score as metrics.

- In binary condition evaluation, binary f1-score is applied. Since F1 score evaluates with regard to the positive categories, you may want to make the 'case' condition labeled as positive. For doing this, you could simply make in the -smp-info file with the 'control' condition appears in the first row, so that the control group will be labelled as 0, and the case group will be labeled as 1.
- In multi-condition evaluation, micro f1-score is applied. In this case, all classes are evaluated equally (no preference of case over control). Furtherly, because each sample is supposed to have only one label, the micro f1-score here is actually equivalent to the accuracy. (cf. [here](https://towardsdatascience.com/multi-class-metrics-made-simple-part-ii-the-f1-score-ebe8b2c2ca1))

For more information about binary/micro f1-score in MLPack, please find [here](https://mlpack.org/doc/mlpack-3.1.0/doxygen/classmlpack_1_1cv_1_1F1.html).

For more information about category-specific metrics, please find [here](https://github.com/jmgirard/mReliability/wiki/Specific-agreement-coefficient).

### Input Count Matrix for KaMRaT

The input count matrix should be in .tsv format, where values are separated by tabulation character.

In the matrix, features are presented as rows, and samples as columns. The first column in matrix should always be the feature column. 

KaMRaT accepts extra columns representing non-count values (e.g. feature's p-value, score, etc.), in this case, please provide a smp-info file indicating which columns are the count columns if necessary (e.g. when merging with intervention, ranking, etc.).

### Output Count Matrix by KaMRaT

The output count matrix is also in .tsv format, where values are separated by tabulation character.

In the matrix, the reduced features are presented as rows, and the columns are in same order as the input.

For KaMRaT rank, an extra column named 'score' is inserted as the second column.

## KaMRaT Execution

You can run KaMRaT under the pattern:

```bash
/..path_to_KaMRaT../bin/kamrat <CMD> [option] input_table > output_table # <CMD> can be one from filter, mask, merge, norm, rank
```

If using Singularity, please add a preceding command ```singularity exec``` before the KaMRaT command line, and you can ignore the leading path to KaMRaT ```bin/``` folder, i.e.

```bash
singularity exec kamrat <CMD> [option] input_table > output_table # <CMD> can be one from filter, mask, merge, norm, rank
```

In order to avoid being wordy, the following guidance only shows KaMRaT execution without Singularity. In the case of using Singularity, please add the preceding command and ignore the leading path to kamrat executable file, as mentioned above.

### KaMRaT Modules

KaMRaT contains 5 modules:

- Filter takes a k-mer count matrix as input, then output the submatrix of k-mers that satisfy the given expressed/silent criteria
- Mask takes a k-mer count matrix and a contig fasta file as input, then output the submatrix of k-mers that appear in contig sequences
- Merge takes a k-mer count matrix as input, then output a contig count matrix by merging k-mers according to their overlap (with/without intervention by count vector)
- Norm takes a feature* count matrix as input, then output the normalized feature count matrix
- Rank takes a feature* count matrix as input, scores and ranks the features on their association with sample labels, and output scores with counts for all or top-N features

\* The feature count matrix can be not only k-mer count matrix, but any kind of count matrix once the first column represents the feature (gene, transcript, etc.)

### KaMRaT Helper

KaMRaT's top-level helper is accessible by typing one of these commands:

```bash
/..path_to_KaMRaT../bin/kamrat
/..path_to_KaMRaT../bin/kamrat -h
/..path_to_KaMRaT../bin/kamrat -help
```

Helpers of each KaMRaT modules are accessible via:

```bash
# <CMD> can be one from filter, mask, merge, norm, rank #
/..path_to_KaMRaT../bin/kamrat <CMD> -h
/..path_to_KaMRaT../bin/kamrat <CMD> -help
```

### KaMRaT filter

Filter usage:

```text
/..path_to_KaMRaT../bin/kamrat filter -express-name STR -silent-name STR -express-thres INT_REC:INT_ABD -silent-thres INT_REC:INT_ABD [-smp-info STR] KMER_COUNT_TAB_PATH > KMER_COUNT_SUB_TAB_PATH
```

Filter parameters:

```text
-h,-help                      Print the helper
-express-name STR             Expressed samples as filter [mandatory]
                                  all         for considering all samples
                                  rest        for considering samples not related with -silent-name
                                  STR:cond    for considering samples in the condition indicated by STR
                                  STR:smp     for considering the only sample indicated by STR
-silent-name STR              Silent samples as filter [mandatory]
                                  all         for considering all samples
                                  rest        for considering samples not related with -express-name
                                  STR:cond    for considering samples in the condition indicated by STR
                                  STR:smp     for considering the only sample indicated by STR
-express-thres INT_R:INT_A    Expressed filter threshold [mandatory]
                                  keep features which have at least (>=) INT_R samples with count >= INT_A (only consider columns related with -express-name)
                                  if -express-name indicates a sample name, INT_R must be 1
-silent-thres INT_R:INT_A     Silent filter threshold [mandatory]
                                  keep features which have at least (>=) INT_R samples with count <= INT_A (only consider columns related with -silent-name)
                                  if -silent-name indicates a sample name, INT_R must be 1
-smp-info STR                 Path to sample-condition file, without header line
                                  if absent, all columns except the first are regarded as samples and labeled as \"all\"
```

### KaMRaT mask

Mask usage:

```text
/..path_to_KaMRaT../bin/kamrat mask -klen INT -fasta STR [-unstrand] [-reverse-mask] KMER_COUNT_TAB_PATH > KMER_COUNT_SUB_TAB_PATH
```

Mask parameters:

```text
-h,-help         Print the helper
-klen INT        k-mer length [mandatory, default value: 31]
-fasta STR       Sequence fasta file as mask [mandatory]
-unstrand        If consider also the reverse complement of k-mers in sequence fasta
-reverse-mask    Reversely mask, i.e. to select (rather than to remove) the k-mers in the sequence fasta file
```

### KaMRaT merge

Merge usage:

```text
/..path_to_KaMRaT../bin/kamrat merge -klen INT [-min-overlap INT] [-unstrand] [-smp-info STR] [-interv-method STR] [-quant STR] [-rep-name STR] [-disk] [-idx-dir STR] KMER_COUNT_TAB_PATH > CONTIG_COUNT_TAB_PATH
```

Merge parameters:

```text
-h,-help              Print the helper
-klen INT             k-mer length [mandatory, max value: 32]
-unstrand             If the k-mers are generated from non-stranded RNA-seq data
-min-overlap INT      Min assembly overlap [max value: k, default value: floor(k/2)]
-smp-info STR         Path to sample-condition file, without header line
                          if absent, all columns except the first are regarded as samples
-interv-method STR    Intervention method [chosen from {none, pearson, spearman, mac}, default value: none]
                          threshold can be precised after a ':' symbol [default: pearson:0.61, spearman:0.56, mac:0.25]
-quant STR            Quantification mode [value chosen from {rep, mean}, defalut value: rep]
-rep-name STR         Column name for representative value
                          if absent, input order of k-mers are taken for choosing representative k-mer for each contig
-disk                 Query on disk [default value: false]
-idx-dir STR          Count index directory path [default value: ./]
```

Intervention method:

```text
none:        No intervention, merge the pre-contigs once they have unique overlap
mac:         Mean absolute contrast between pre-contigs to merge, mac(c1, c2) = mean(abs(c1-c2)./(c1+c2)), where c1, c2 are count vectors of two pre-contigs and './' is division for each components.
pearson:     Pearson correlation between pre-condtigs to merge
spearman:    Spearman correlation between pre-contigs to merge
```

Count vector for output:

```text
rep:     representative* k-mer's sample count vector and other non-count values is output for each contig
mean:    average sample count vector of all composite k-mers and other non-count values of the representative* k-mer is output for each contig

* if no -rep-name is provided, the first input k-mer will be taken as the representative k-mer for each contig
```

### KaMRaT norm

Norm usage:

```text
/..path_to_KaMRaT../bin/kamrat norm -base CHAR [-smp-info STR] [-ln] [-smp-sum STR] COUNT_TAB_PATH > NORM_COUNT_TAB_PATH
```

Norm parameters:

```text
-h,-help         Print the helper
-base CHAR       Base for normalization [MANDATORY]
                     B: count per billion
                     M: count per million
                     K: count per thousand
-smp-info STR    Path to sample-condition file, without header line
                     if absent, all columns except the first are regarded as samples and will be normalized
-ln              If apply ln(x + 1) transformation after normalization
                     hint: please remember to unlog for original counts after kamratReduce or other analysis
-smp-sum STR     Path for outputing sample sum [default value: ./sample_sum.tsv]
```

### KaMRaT rank

Rank usage:

```text
/..path_to_KaMRaT../bin/kamrat rank [-smp-info STR] [-eval-method STR] [-sort-mode STR] [-top-num INT] [-ln] COUNT_TAB_PATH > COUNT_SUB_TAB_PATH
```

Rank parameters:

```text
-h,-help             Print the helper
-smp-info STR        Path to sample-condition file, without header line
                          if absent, all columns except the first are regarded as samples
-score-method STR    Evaluation method to use and its parameter, seperated by \':\' [default value: sd]
                         sd            Standard deviation
                         rsd           Relative standard deviation
                         ttest         T-test adjusted p-value between conditions
                         es            Effect size between conditions
                         lfc:mean      Log2 fold change by group mean, 'mean' can be omitted by default
                         lfc:median    Log2 fold change by group median
                         nb:n_fold     Naive Bayes classification [default n_fold = 1 (no cross-validation)]
                                           if n_fold = 0, leave-one-out cross-validation is applied
                                           if n_fold = 1, no cross-validation is applied, features are evaluated by training and testing on the whole datset
                                           if n_fold >= 2, n-fold cross-validation is applied
                         rg:n_fold     Classification by regression (logistic regression for binary conditions, softmax regression for multiple conditions) [default n_fold = 1]
                                           if n_fold = 0, leave-one-out cross-validation is applied
                                           if n_fold = 1, no cross-validation is applied, features are evaluated by training and testing on the whole datset
                                           if n_fold >= 2, n-fold cross-validation is applied
                         user:NAME     User-defined method, NAME indicates the score column in the k-mer count table
-sort-mode STR       Mode for sorting features by score, default value depends on evaluation method
                         dec        sort by decreasing order                              [as default value for sd, rsd, nb, rg, user:name]
                         dec:abs    sort by decreasing order but on the absolute value    [as default value for es, lfc:mean, lfc:median]
                         inc        sort by increasing order                              [as default value for ttest]
                         inc:abs    sort by increasing order but on the absolute value
-top-num INT         Number of top features to output
-ln                  Apply ln(x + 1) transformation BEFORE score estimation [default value: false]
                         note: this applies ONLY for score estimation, will NOT affect output counts
```

## Software/Library Citations

Armadillo:

- Conrad Sanderson and Ryan Curtin. Armadillo: a template-based C++ library for linear algebra. Journal of Open Source Software, Vol. 1, pp. 26, 2016.
- Conrad Sanderson and Ryan Curtin. A User-Friendly Hybrid Sparse Matrix Class in C++. Lecture Notes in Computer Science (LNCS), Vol. 10931, pp. 422-430, 2018.

[Boost C++ Library](https://www.boost.org/)

DE-kupl: Audoux, J., Philippe, N., Chikhi, R. et al. DE-kupl: exhaustive capture of biological variation in RNA-seq data through k-mer decomposition. Genome Biol 18, 243 (2017).

MLPack: R.R. Curtin, M. Edel, M. Lozhnikov, Y. Mentekidis, S. Ghaisas, S. Zhang. mlpack 3: a fast, flexible machine learning library. Journal of Open Source Software 3:26, 2018.
