#ifndef KMEREVALUATE_UTILS_H
#define KMEREVALUATE_UTILS_H

#include <iostream>
#include <string>

#define CLASS_NUM 2
#define METHOD_NB "nb"
#define METHOD_LR "lr"
#define METHOD_SD "sd"
#define METHOD_SDC "sdc"
#define METHOD_RSD "rsd"
#define METHOD_TTEST "ttest"
#define METHOD_ES "es"
#define METHOD_USER "user"
#define SORT_DEC "dec"
#define SORT_INC "inc"
#define SORT_ABS "abs"

typedef double (*evalfuncptr_t)(const arma::mat &, const arma::Row<size_t> &, size_t, size_t);

inline void PrintHelper()
{
    std::cerr << "=====> kamratReduce Helper <=====" << std::endl;
    std::cerr << "[Usage]         kamratReduce [-d samp_info_path -m eval_method:fold_num -s sort_mode -n top_num] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Parameters]    -d STRING    Path to sample-condition or sample file, without header line" << std::endl
              << "                             if absent, all except the first column in k-mer count table will be regarded as samples" << std::endl;
    std::cerr << "                -m STRING    Evaluation method name and fold number for cross-validation (if needed), seperated by \':\'" << std::endl;
    std::cerr << "                -s STRING    Sorting mode, default value depends on evaluation method (c.f [SORT MODE])" << std::endl;
    std::cerr << "                -n INT       Number of top features to print" << std::endl
              << std::endl;
    std::cerr << "[Evaluation]    nb           Naive Bayes classification between conditions" << std::endl;
    std::cerr << "                lr           Logistic regression (slower than Naive Bayes) between conditions" << std::endl;
    std::cerr << "                sd           Standard derivation" << std::endl;
    std::cerr << "                rsd          Relative standard derivation" << std::endl;
    std::cerr << "                sdc          Standard deviation contrast between conditions" << std::endl;
    std::cerr << "                ttest        T-test between conditions" << std::endl;
    std::cerr << "                es           Effect size between conditions" << std::endl;
    std::cerr << "                user:name    User-defined method, where name indicates a column in the k-mer count table" << std::endl
              << std::endl;
    std::cerr << "[SORT MODE]     dec          Sorting by decreasing order                              (as default for nb, lr, sd, rsd, user:name)" << std::endl;
    std::cerr << "                dec:abs      Sorting by decreasing order but on the absolute value    (as default for sdc, es)" << std::endl;
    std::cerr << "                inc          Sorting by increasing order                              (as default for ttest)" << std::endl;
    std::cerr << "                inc:abs      Sorting by increasing order but on the absolute value" << std::endl
              << std::endl;
}

inline void PrintParameterInfo(const std::string &kmer_count_path,
                               const std::string &sample_info_path,
                               const std::string &eval_method,
                               const size_t fold_num,
                               const std::string &sort_mode,
                               const size_t top_num)
{
    std::cerr << "=====> PARAMETER INFO <=====" << std::endl;
    std::cerr << "\tk-mer count path: " << kmer_count_path << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "\tsample info path: " << sample_info_path << std::endl;
    }
    std::cerr << "\tevaluation method: " << eval_method << std::endl;
    if (eval_method == METHOD_NB || eval_method == METHOD_LR)
    {
        std::cerr << "\tk-fold number: " << fold_num << std::endl;
    }
    if (!sort_mode.empty())
    {
        std::cerr << "\tsorting mode: " << sort_mode << std::endl;
    }
    if (top_num > 0)
    {
        std::cerr << "\tnumber of k-mers to select: " << top_num << std::endl;
    }
    std::cerr << std::endl;
    if (eval_method == METHOD_RSD)
    {
        std::cerr << "[ATTENTION] Please note that the relative standard deviation is only suitable for non-normalized counts" << std::endl;
    }
}

inline void CheckInputError(const bool arg_finished)
{
    if (arg_finished)
    {
        std::cerr << "ERROR: k-mer count table path is mandatory." << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

#endif //KMEREVALUATE_UTILS_H