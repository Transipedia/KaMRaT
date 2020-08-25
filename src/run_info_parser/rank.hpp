#ifndef KAMRAT_RUNINFOPARSER_RANK_HPP
#define KAMRAT_RUNINFOPARSER_RANK_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintRankHelper()
{
    std::cerr << "[USAGE]    kamrat rank [-smp-info STR] [-eval-method STR] [-sort-mode STR] [-top-num INT] [-ln] COUNT_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]        -smp-info STR        Path to sample-condition or sample file, without header line" << std::endl
              << "                                         if absent, all columns except the first in the count table are regarded as sample" << std::endl;
    std::cerr << "                -score-method STR    Evaluation method to use and its parameter, seperated by \':\' (cf. [EVAL. METHOD])" << std::endl;
    std::cerr << "                -sort-mode STR       Mode for score sorting, default value depends on evaluation method (cf. [SORT MODE])" << std::endl;
    std::cerr << "                -top-num INT         Number of top features to select" << std::endl;
    std::cerr << "                -ln                  Apply ln(x + 1) transformation BEFORE score estimation [false]" << std::endl
              << "                                         this applies only for score estimation, will NOT affect output counts" << std::endl
              << std::endl;
    std::cerr << "[EVAL. METHOD]  sd                  Standard deviation (default method)" << std::endl;
    std::cerr << "                rsd                 Relative standard deviation" << std::endl;
    std::cerr << "                ttest               T-test between conditions" << std::endl;
    std::cerr << "                es                  Effect size between conditions" << std::endl;
    std::cerr << "                lfc:mean            Log2 fold change by group mean, 'mean' can be omitted by default" << std::endl;
    std::cerr << "                lfc:median          Log2 fold change by group median" << std::endl;
    std::cerr << "                nb:n_fold           Naive Bayes classification, default n_fold = 2" << std::endl;
    std::cerr << "                rg:n_fold           Classification by regression, default n_fold = 2" << std::endl
              << "                                        if nb_condition = 2, logistic regression is used" << std::endl
              << "                                        if nb_condition > 2, softmax regression is used" << std::endl;
    std::cerr << "                user:name           User-defined method, where name indicates a column in the k-mer count table" << std::endl
              << std::endl;
    std::cerr << "[SORT MODE]     dec                 Sorting by decreasing order                              (as default for sd, rsd, nb, lr, user:name)" << std::endl;
    std::cerr << "                dec:abs             Sorting by decreasing order but on the absolute value    (as default for es, lfc:mean, lfc:median)" << std::endl;
    std::cerr << "                inc                 Sorting by increasing order                              (as default for ttest)" << std::endl;
    std::cerr << "                inc:abs             Sorting by increasing order but on the absolute value" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &kmer_count_path,
                         const std::string &sample_info_path,
                         const std::string &score_method,
                         const std::string &score_cmd,
                         const std::string &sort_mode,
                         const size_t nb_fold,
                         const size_t nb_sel,
                         const bool ln_transf)
{
    std::cerr << "k-mer count path:                             " << kmer_count_path << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample info path:                             " << sample_info_path << std::endl;
    }
    std::cerr << "Evaluation method:                            " << score_method << std::endl;
    if (nb_fold >= 2)
    {
        std::cerr << "    Fold number = " << nb_fold << std::endl;
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
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &sample_info_path,
                         std::string &score_method,
                         std::string &score_cmd,
                         std::string &sort_mode,
                         size_t &nb_sel,
                         bool &ln_transf,
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
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
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
        else
        {
            PrintRankHelper();
            throw std::domain_error("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintRankHelper();
        throw std::domain_error("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];
    if (!sort_mode.empty() && SORT_MODE_UNIV.find(sort_mode) == SORT_MODE_UNIV.cend())
    {
        PrintRankHelper();
        throw std::domain_error("unknown sort mode: " + sort_mode);
    }
}

#endif //KAMRAT_RUNINFOPARSER_RANK_HPP
