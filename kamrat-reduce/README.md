# kamratReduce Helper

## Usage

    kamratReduce [-d samp_info_path -m eval_method:fold_num -s sort_mode -n top_num] kmer_count_path

## Parameters

    -d STRING    Path to sample-condition or sample file, without header line
                 if absent, all except the first column in k-mer count table will be regarded as samples
    -m STRING    Evaluation method name and fold number for cross-validation (if needed), seperated by ':'
    -s STRING    Sorting mode, default value depends on evaluation method (c.f [SORT MODE])
    -n INT       Number of top features to print

## Evaluation Methods

    nb           Naive Bayes classification between conditions
    lr           Logistic regression (slower than Naive Bayes) between conditions
    sd           Standard deviation
    rsd          Relative standard deviation
    sdc          Standard deviation contrast between conditions
    ttest        T-test between conditions
    user:name    User-defined method, where name indicates a column in the k-mer count table

## Sorting Modes

    dec          Sorting by decreasing order                              (as default for nb, lr, sd, rsd, user:name)
    dec:abs      Sorting by decreasing order but on the absolute value    (as default for sdc)
    inc          Sorting by increasing order                              (as default for ttest)
    inc:abs      Sorting by increasing order but on the absolute value
