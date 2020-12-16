#ifndef KAMRAT_RUNINFOPARSER_RANK_HPP
#define KAMRAT_RUNINFOPARSER_RANK_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintRankHelper()
{
    std::cerr << "[USAGE]   kamrat rank -idx-path STR -nf-path STR [-options] FEATURE_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]        -h,-help             Print the helper " << std::endl;
    std::cerr << "                -idx-path STR        Temporary file path for saving count index, mandatory" << std::endl;
    std::cerr << "                -nf-path             Output path for nomalization factor, mandatory" << std::endl;
    std::cerr << "                -smp-info STR        Path to sample-condition or sample file, without header line" << std::endl
              << "                                         if absent, all columns except the first in the count table are regarded as sample" << std::endl;
    std::cerr << "                -score-method STR    Evaluation method to use and its parameter, seperated by \':\' (cf. [EVAL. METHOD])" << std::endl;
    std::cerr << "                -sort-mode STR       Mode for score sorting, default value depends on evaluation method (cf. [SORT MODE])" << std::endl;
    std::cerr << "                -top-num INT         Number of top features to select" << std::endl;
    std::cerr << "                -ln                  Apply ln(x + 1) transformation BEFORE score estimation [false]" << std::endl;
    std::cerr << "                -standardize         Standarize count vector BEFORE score estimation [false]" << std::endl;
    std::cerr << "                -out-path STR        Output table path [default: output to screen]" << std::endl
              << "                                         the output counts are same as the input counts," << std::endl
              << "                                         normalization, log transformation, or standardization affects only score evaluation, not output counts" << std::endl
              << std::endl;
    std::cerr << "[EVAL. METHOD]  sd                   Standard deviation (default method)" << std::endl;
    std::cerr << "                rsd                  Relative standard deviation" << std::endl;
    std::cerr << "                ttest                T-test between conditions (ln transformation is required)" << std::endl;
    std::cerr << "                snr                  Signal-to-noise ratio between conditions" << std::endl;
    std::cerr << "                lfc:mean             Log2 fold change by group mean, 'mean' can be omitted by default" << std::endl;
    std::cerr << "                lfc:median           Log2 fold change by group median" << std::endl;
    std::cerr << "                nb:n_fold            F1-score with naive Bayes classification [default n_fold = 1]" << std::endl
              << "                                         if n_fold = 0, leave-one-out cross-validation is applied" << std::endl
              << "                                         if n_fold = 1, no cross-validation is applied, features are evaluated by training and testing on the whole datset" << std::endl
              << "                                         if n_fold >= 2, n-fold cross-validation is applied" << std::endl;
    std::cerr << "                rg:n_fold            F1-score with regression classification [default n_fold = 1]" << std::endl
              << "                                         if n_fold = 0, leave-one-out cross-validation is applied" << std::endl
              << "                                         if n_fold = 1, no cross-validation is applied, features are evaluated by training and testing on the whole datset" << std::endl
              << "                                         if n_fold >= 2, n-fold cross-validation is applied" << std::endl;
    std::cerr << "                svm                  F1-score with SVM classification (standardization is required)" << std::endl;
    std::cerr << "                user:name            User-defined method, where name indicates a column in the k-mer count table" << std::endl
              << std::endl;
    std::cerr << "[SORT MODE]     dec                  Sorting by decreasing order                              (as default for sd, rsd, nb, lr, user:name)" << std::endl;
    std::cerr << "                dec:abs              Sorting by decreasing order but on the absolute value    (as default for es, lfc:mean, lfc:median)" << std::endl;
    std::cerr << "                inc                  Sorting by increasing order                              (as default for ttest)" << std::endl;
    std::cerr << "                inc:abs              Sorting by increasing order but on the absolute value" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &kmer_count_path,
                         const std::string &idx_path,
                         const std::string &nf_path,
                         const std::string &smp_info_path,
                         const std::string &score_method,
                         const std::string &score_cmd,
                         const std::string &sort_mode,
                         const size_t nb_fold,
                         const size_t nb_sel,
                         const bool ln_transf,
                         const bool standardize,
                         const std::string &out_path)
{
    std::cerr << "k-mer count path:                             " << kmer_count_path << std::endl;
    std::cerr << "k-mer count index path:                       " << idx_path << std::endl;
    if (!smp_info_path.empty())
    {
        std::cerr << "Sample info path:                             " << smp_info_path << std::endl;
    }
    std::cerr << "Evaluation method:                            " << score_method << std::endl;
    if (score_method == "naivebayes.f1" || score_method == "regression.f1")
    {
        if (nb_fold == 0)
        {
            std::cerr << "    Leave-one-out cross-validation";
        }
        else if (nb_fold == 1)
        {
            std::cerr << "    No cross-validation, train and test on the whole dataset";
        }
        else
        {
            std::cerr << "    Cross-validation with fold number = " << nb_fold;
        }
        std::cerr << std::endl;
    }
    else if (score_method == "log2fc")
    {
        std::cerr << "    Value for comparison = " << score_cmd << " of group" << std::endl;
    }
    else if (score_method == "user")
    {
        std::cerr << "    Score column name = " << score_cmd << std::endl;
    }
    if (!sort_mode.empty())
    {
        std::cerr << "Sorting mode:                                 " << sort_mode << std::endl;
    }
    std::cerr << "Number of feature to output (0 for all):      " << nb_sel << std::endl;
    std::cerr << "Ln(x + 1) for score estiamtion:               " << (ln_transf ? "On" : "Off") << std::endl;
    std::cerr << "Standardize for score estimation:             " << (standardize ? "On" : "Off") << std::endl;
    std::cerr << "Nomalization factor to path:                  " << nf_path << std::endl;
    if (!out_path.empty())
    {
        std::cerr << "Output path:                                  " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &idx_path,
                         std::string &nf_path,
                         std::string &smp_info_path,
                         std::string &score_method,
                         std::string &score_cmd,
                         std::string &sort_mode,
                         size_t &nb_sel,
                         bool &ln_transf,
                         bool &standardize,
                         std::string &out_path,
                         std::string &kmer_count_path)
{
    int i_opt(1);
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h" || arg == "-help")
        {
            PrintRankHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-idx-path" && i_opt + 1 < argc)
        {
            idx_path = argv[++i_opt];
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            smp_info_path = argv[++i_opt];
        }
        else if (arg == "-score-method" && i_opt + 1 < argc)
        {
            score_method = argv[++i_opt];
            SubCommandParser(score_method, score_cmd);
        }
        else if (arg == "-sort-mode" && i_opt + 1 < argc)
        {
            sort_mode = argv[++i_opt];
        }
        else if (arg == "-top-num" && i_opt + 1 < argc)
        {
            nb_sel = atoi(argv[++i_opt]);
        }
        else if (arg == "-ln")
        {
            ln_transf = true;
        }
        else if (arg == "-standardize")
        {
            standardize = true;
        }
        else if (arg == "-nf-path" && i_opt + 1 < argc)
        {
            nf_path = argv[++i_opt];
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else
        {
            PrintRankHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintRankHelper();
        throw std::invalid_argument("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];
    if (idx_path.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("temporary index file path is mandatory");
    }
    if (!sort_mode.empty() && SORT_MODE_UNIV.find(sort_mode) == SORT_MODE_UNIV.cend())
    {
        PrintRankHelper();
        throw std::invalid_argument("unknown sort mode: " + sort_mode);
    }
    if (nf_path.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("Path for normalization factor is mandatory");
    }
    if ((score_method == "sd" || score_method == "rsd") && standardize)
    {
        throw std::invalid_argument("Standard deviation score should not be applied on standardized counts\n"
                                    "And the forced running is not permitted either");
    }
    if (score_method == "svm" && !standardize && score_cmd != "F")
    {
        throw std::invalid_argument("Standardization is required for SVM classification\n"
                                    "Put -score-method svm:F to force to run");
    }
    if (score_method == "ttest" && standardize)
    {
        throw std::invalid_argument("Ttest is not compatible with standardized counts: log applied on negative values\n"
                                    "And the forced running is not permitted either");
    }
    if (score_method == "ttest" && !ln_transf && score_cmd != "F")
    {
        throw std::invalid_argument("Ttest requires ln(x + 1) transformation\n"
                                    "Put -score-method ttest:F to force to run");
    }
}

#endif //KAMRAT_RUNINFOPARSER_RANK_HPP
