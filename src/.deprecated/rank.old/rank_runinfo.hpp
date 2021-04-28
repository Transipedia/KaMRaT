#ifndef KAMRAT_RANK_RANKRUNINFO_HPP
#define KAMRAT_RANK_RANKRUNINFO_HPP

#include <iostream>
#include <string>
#include <memory>

#include "scorer.hpp"

const void ParseScorer(std::unique_ptr<Scorer> &scorer, std::string &score_method, const bool to_ln, const bool to_standardize)
{
    std::string score_cmd;
    size_t split_pos = score_method.find(":");
    if (split_pos != std::string::npos)
    {
        score_cmd = score_method.substr(split_pos + 1);
        score_method = score_method.substr(0, split_pos);
    }
    if (score_method == "rsd")
    {
        if (to_standardize)
        {
            throw std::invalid_argument("Standard deviation score should not be applied on standardized counts\n"
                                        "And the forced running is not permitted either (divided by 0)");
        }
        scorer = std::make_unique<Scorer>(ScoreMethodCode::kRelatSD);
    }
    else if (score_method == "ttest")
    {
        if (!to_ln && score_cmd != "F")
        {
            throw std::invalid_argument("Ttest requires ln(x + 1) transformation, put -score-method ttest:F to force to run");
        }
        scorer = std::make_unique<Scorer>(ScoreMethodCode::kTtest);
    }
    else if (score_method == "snr")
    {
        scorer = std::make_unique<Scorer>(ScoreMethodCode::kSNR);
    }
    else if (score_method == "lrc")
    {
        if (!score_cmd.empty())
        {
            scorer = std::make_unique<Scorer>(ScoreMethodCode::kLogitReg, std::stoul(score_cmd)); // C++: typedef unsigned long size_t
        }
        else
        {
            scorer = std::make_unique<Scorer>(ScoreMethodCode::kLogitReg, 1); // by default, evaluate features without cross-validation
        }
    }
    else if (score_method == "nbc")
    {
        if (!score_cmd.empty())
        {
            scorer = std::make_unique<Scorer>(ScoreMethodCode::kNaiveBayes, std::stoul(score_cmd)); // C++: typedef unsigned long size_t
        }
        else
        {
            scorer = std::make_unique<Scorer>(ScoreMethodCode::kNaiveBayes, 1); // by default, evaluate features without cross-validation
        }
    }
    else if (score_method == "svm")
    {
        if (!to_standardize && score_cmd != "F")
        {
            throw std::invalid_argument("Standardization is required for SVM classification, put -score-method svm:F to force to run");
        }
        scorer = std::make_unique<Scorer>(ScoreMethodCode::kSVM);
    }
    else
    {
        if (score_cmd == "dec")
        {
            scorer = std::make_unique<Scorer>(score_method, SortModeCode::kDec);
        }
        else if (score_cmd == "decabs")
        {
            scorer = std::make_unique<Scorer>(score_method, SortModeCode::kDecAbs);
        }
        else if (score_cmd == "inc")
        {
            scorer = std::make_unique<Scorer>(score_method, SortModeCode::kInc);
        }
        else if (score_cmd == "incabs")
        {
            scorer = std::make_unique<Scorer>(score_method, SortModeCode::kIncAbs);
        }
        else if (score_cmd.empty())
        {
            throw std::invalid_argument("user-defined ranking with column (" + score_method + ") but without indicating sorting mode");
        }
        else
        {
            throw std::invalid_argument("unknown sorting mode: " + score_cmd);
        }
        
    }
}

inline void PrintRankHelper()
{
    std::cerr << "[USAGE]   kamrat rank -idx-path STR -nf-path STR [-options] FEATURE_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]        -h,-help             Print the helper " << std::endl;
    std::cerr << "                -idx-path STR        Temporary file path for saving count index [MANDATORY]" << std::endl;
    std::cerr << "                -nf-path             Output path for nomalization factor [MANDATORY]" << std::endl;
    std::cerr << "                -smp-info STR        Path to sample-condition or sample file, without header line" << std::endl
              << "                                         if absent, all columns except the first in the count table are regarded as sample" << std::endl;
    std::cerr << "                -score-method STR    Evaluation method to use and its parameter, seperated by \':\' (cf. [EVAL. METHOD])" << std::endl;
    std::cerr << "                -top-num INT         Number of top features to select" << std::endl;
    std::cerr << "                -ln                  Apply ln(x + 1) transformation for score estimation [false]" << std::endl;
    std::cerr << "                -standardize         Standarize count vector for score estimation [false]" << std::endl;
    std::cerr << "                -no-norm             Estimate scores with raw count, do NOT apply normalization" << std::endl;
    std::cerr << "                -out-path STR        Output table path [default: output to screen]" << std::endl
              << "                                         the output counts are same as the input counts," << std::endl
              << "                                         normalization, log transformation, and standardization affect score evaluation, but not output counts" << std::endl
              << std::endl;
    std::cerr << "[EVAL. METHOD]  rsd                  Relative standard deviation" << std::endl;
    std::cerr << "                ttest                T-test between conditions (ln transformation is required)" << std::endl;
    std::cerr << "                snr                  Signal-to-noise ratio between conditions" << std::endl;
    std::cerr << "                lrc:n_fold           F1-score with regression classification [default n_fold = 1]" << std::endl
              << "                                         if n_fold = 0, leave-one-out cross-validation is applied" << std::endl
              << "                                         if n_fold = 1, evaluation without cross-validation, training and testing on the whole datset" << std::endl
              << "                                         if n_fold >= 2, n-fold cross-validation is applied" << std::endl;
    std::cerr << "                nbc:n_fold           F1-score with naive Bayes classification [default n_fold = 1]" << std::endl
              << "                                         if n_fold = 0, leave-one-out cross-validation is applied" << std::endl
              << "                                         if n_fold = 1, evaluation without cross-validation, training and testing on the whole datset" << std::endl
              << "                                         if n_fold >= 2, n-fold cross-validation is applied" << std::endl;
    std::cerr << "                svm                  Hinge-loss function on SVM classification (standardization is required)" << std::endl;
    std::cerr << "                colname:sort_mode    User-defined method, where name indicates a column in the k-mer count table" << std::endl
              << "                                         sore_mode can be:    dec       Sorting by decreasing order" << std::endl
              << "                                                              decabs    Sorting by decreasing order but on the absolute value" << std::endl
              << "                                                              inc       Sorting by increasing order" << std::endl
              << "                                                              incabs    Sorting by increasing order but on the absolute value" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &kmer_count_path,
                         const std::string &idx_path,
                         const std::string &nf_path,
                         const std::string &smp_info_path,
                         const std::unique_ptr<Scorer> &scorer,
                         const size_t nb_sel,
                         const bool ln_transf,
                         const bool standardize,
                         const bool no_norm,
                         const std::string &out_path)
{
    std::cerr << "k-mer count path:                             " << kmer_count_path << std::endl;
    std::cerr << "k-mer count index path:                       " << idx_path << std::endl;
    std::cerr << "Nomalization factor to path:                  " << nf_path << std::endl;
    if (!smp_info_path.empty())
    {
        std::cerr << "Sample info path:                             " << smp_info_path << std::endl;
    }
    const ScoreMethodCode score_method_code = scorer->GetScoreMethodCode();
    std::cerr << "Evaluation method:                            " << kScoreMethodName[score_method_code] << std::endl;
    if (score_method_code == ScoreMethodCode::kNaiveBayes || score_method_code == ScoreMethodCode::kLogitReg)
    {
        size_t nb_fold = scorer->GetNbFold();
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
    else if (score_method_code == ScoreMethodCode::kUser)
    {
        std::cerr << "    Score column name = " << scorer->GetRepColname() << std::endl;
    }
    std::cerr << "Sorting mode:                                 " << kSortModeName[scorer->GetSortModeCode()] << std::endl;
    std::cerr << "Number of feature to output (0 for all):      " << nb_sel << std::endl;
    std::cerr << "Ln(x + 1) for score estiamtion:               " << (ln_transf ? "On" : "Off") << std::endl;
    std::cerr << "Standardize for score estimation:             " << (standardize ? "On" : "Off") << std::endl;
    std::cerr << "Feature evaluation:                           on " << (no_norm ? "raw" : "normalized") << " counts" << std::endl;
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
                         std::unique_ptr<Scorer> &scorer,
                         size_t &nb_sel,
                         bool &to_ln,
                         bool &to_standardize,
                         bool &no_norm,
                         std::string &out_path,
                         std::string &kmer_count_path)
{
    std::string score_method;
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
        else if (arg == "-nf-path" && i_opt + 1 < argc)
        {
            nf_path = argv[++i_opt];
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            smp_info_path = argv[++i_opt];
        }
        else if (arg == "-score-method" && i_opt + 1 < argc)
        {
            score_method = argv[++i_opt];
        }
        else if (arg == "-top-num" && i_opt + 1 < argc)
        {
            nb_sel = atoi(argv[++i_opt]);
        }
        else if (arg == "-ln")
        {
            to_ln = true;
        }
        else if (arg == "-standardize")
        {
            to_standardize = true;
        }
        else if (arg == "-no-norm")
        {
            no_norm = true;
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
    if (nf_path.empty() && !no_norm)
    {
        PrintRankHelper();
        throw std::invalid_argument("path for normalization factor is mandatory unless with -no-norm");
    }
    else if (!nf_path.empty() && no_norm)
    {
        PrintRankHelper();
        throw std::invalid_argument("no need for giving normalization factor path with -no-norm");
    }

    ParseScorer(scorer, score_method, to_ln, to_standardize);
}

#endif //KAMRAT_RANK_RANKRUNINFO_HPP
