#ifndef KAMRAT_RANK_RANKRUNINFO_HPP
#define KAMRAT_RANK_RANKRUNINFO_HPP

#include <iostream>
#include <string>
#include <memory>

// const void ParseScorer(std::unique_ptr<Scorer> &scorer, std::string &score_method, const bool ln_transf, const bool standardize)
// {
//     std::string score_cmd;
//     size_t split_pos = score_method.find(":");
//     if (split_pos != std::string::npos)
//     {
//         score_cmd = score_method.substr(split_pos + 1);
//         score_method = score_method.substr(0, split_pos);
//     }
//     if (score_method == "rsd")
//     {
//         if (standardize)
//         {
//             throw std::invalid_argument("Standard deviation score should not be applied on standardized counts\n"
//                                         "And the forced running is not permitted either (divided by 0)");
//         }
//         scorer = std::make_unique<Scorer>(ScoreMethodCode::kRelatSD);
//     }
//     else if (score_method == "ttest")
//     {
//         if (!ln_transf && score_cmd != "F")
//         {
//             throw std::invalid_argument("Ttest requires ln(x + 1) transformation, put -score-method ttest:F to force to run");
//         }
//         scorer = std::make_unique<Scorer>(ScoreMethodCode::kTtest);
//     }
//     else if (score_method == "snr")
//     {
//         scorer = std::make_unique<Scorer>(ScoreMethodCode::kSNR);
//     }
//     else if (score_method == "lrc")
//     {
//         if (!score_cmd.empty())
//         {
//             scorer = std::make_unique<Scorer>(ScoreMethodCode::kLogitReg, std::stoul(score_cmd)); // C++: typedef unsigned long size_t
//         }
//         else
//         {
//             scorer = std::make_unique<Scorer>(ScoreMethodCode::kLogitReg, 1); // by default, evaluate features without cross-validation
//         }
//     }
//     else if (score_method == "nbc")
//     {
//         if (!score_cmd.empty())
//         {
//             scorer = std::make_unique<Scorer>(ScoreMethodCode::kNaiveBayes, std::stoul(score_cmd)); // C++: typedef unsigned long size_t
//         }
//         else
//         {
//             scorer = std::make_unique<Scorer>(ScoreMethodCode::kNaiveBayes, 1); // by default, evaluate features without cross-validation
//         }
//     }
//     else if (score_method == "svm")
//     {
//         if (!standardize && score_cmd != "F")
//         {
//             throw std::invalid_argument("Standardization is required for SVM classification, put -score-method svm:F to force to run");
//         }
//         scorer = std::make_unique<Scorer>(ScoreMethodCode::kSVM);
//     }
//     else
//     {
//         if (score_cmd == "dec")
//         {
//             scorer = std::make_unique<Scorer>(score_method, SortModeCode::kDec);
//         }
//         else if (score_cmd == "decabs")
//         {
//             scorer = std::make_unique<Scorer>(score_method, SortModeCode::kDecAbs);
//         }
//         else if (score_cmd == "inc")
//         {
//             scorer = std::make_unique<Scorer>(score_method, SortModeCode::kInc);
//         }
//         else if (score_cmd == "incabs")
//         {
//             scorer = std::make_unique<Scorer>(score_method, SortModeCode::kIncAbs);
//         }
//         else if (score_cmd.empty())
//         {
//             throw std::invalid_argument("user-defined ranking with column (" + score_method + ") but without indicating sorting mode");
//         }
//         else
//         {
//             throw std::invalid_argument("unknown sorting mode: " + score_cmd);
//         }
//     }
// }

void RankWelcome()
{
    std::cerr << "KaMRaT rank: univariate feature ranking" << std::endl
              << "------------------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintRankHelper()
{
    std::cerr << "[USAGE]    kamrat rank -idx-dir STR -count-mode STR -rank-by STR [-options] FEATURE_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help             Print the helper " << std::endl;
    std::cerr << "            -idxdir STR              Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -rankby STR          Ranking method, mandatory, can be one of: " << std::endl
              << "                                     ttest        adjusted p-value of t-test between conditions (require -ln)" << std::endl
              << "                                     snr          signal-to-noise ratio between conditions" << std::endl
              << "                                     lr:nfold     accuracy by logistic regression classifier" << std::endl
              << "                                     nbc:nfold    accuracy by naive Bayes classifier" << std::endl
              << "                                     svm:nfold    accuracy on SVM classifier" << std::endl;
    std::cerr << "            -with STR1[:STR2]    File indicating features to rank (STR1) and counting mode (STR2)" << std::endl
              << "                                     if not provided, all indexed features are used for ranking" << std::endl
              << "                                     STR2 can be one of [rep, mean, median]" << std::endl;
    std::cerr << "            -design STR          File indicating sample-condition design, without header line" << std::endl
              << "                                     if not provided, all samples are assigned by the same condition" << std::endl
              << "                                     if provided, each row can be either: " << std::endl
              << "                                         sample name, sample condition" << std::endl
              << "                                         sample name, sample condition, sample batch (only for lrc, nbc, and svm)" << std::endl;
    std::cerr << "            -ln                  Apply ln(x + 1) transformation for score estimation [false]" << std::endl;
    std::cerr << "            -standardize         Standarize count vector for score estimation [false]" << std::endl;
    std::cerr << "            -rankonraw           Estimate scores and rank on raw count, without normalization" << std::endl;
    std::cerr << "            -seltop NUM          If NUM > 1, it indicates top number of features to output (treated as integer)" << std::endl
              << "                                 If NUM <= 1, it indicates the ratio of features to output" << std::endl;
    std::cerr << "            -outpath STR         Path of ranking result" << std::endl
              << "                                     if not provided, output to screen" << std::endl;
    std::cerr << "            -withcounts          Output sample count vectors [false]" << std::endl
              << std::endl;
    std::cerr << "[NOTE]      For ranking methods lrc, nbc, and svm, there can be a second univariant cross-validaton option (nfold)" << std::endl
              << "                if nfold = 0, leave-one-out cross-validation" << std::endl
              << "                if nfold = 1, without cross-validation, training and testing on the whole datset" << std::endl
              << "                if nfold > 1, n-fold cross-validation" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir,
                  const std::string &rk_mthd, const size_t nfold,
                  const std::string &with_path, const std::string &count_mode,
                  const std::string &dsgn_path,
                  const bool ln_transf, const bool standardize, const bool no_norm,
                  const float sel_top,
                  const std::string &out_path, const bool with_counts)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                 " << idx_dir << std::endl;
    std::cerr << "Ranking method:               " << rk_mthd;
    if (rk_mthd == "lr" || rk_mthd == "nbc" || rk_mthd == "svm")
    {
        if (nfold == 0)
        {
            std::cerr << ", "
                      << "leave-one-out cross-validation" << std::endl;
        }
        else if (nfold == 1)
        {
            std::cerr << ", "
                      << "no cross-validation, train and test both on all samples" << std::endl;
        }
        else
        {
            std::cerr << ", " << nfold << "-fold cross-validation" << std::endl;
        }
    }
    else
    {
        std::cerr << std::endl;
    }

    std::cerr << "Ranking with:                 " << (with_path.empty() ? "features in index" : "features in " + with_path) << std::endl;
    if (!count_mode.empty())
    {
        std::cerr << "Feature counting mode:        " + count_mode << std::endl;
    }
    if (!dsgn_path.empty())
    {
        std::cerr << "Sample design:                " << dsgn_path << std::endl;
    }
    std::cerr << "ln(x + 1) transformation:     " << (ln_transf ? "On" : "Off") << std::endl;
    std::cerr << "Standardization:              " << (standardize ? "On" : "Off") << std::endl;
    std::cerr << "Rank on raw counts:           " << (no_norm ? "On" : "Off") << std::endl;
    std::cerr << "Selection of top features:    ";
    if (sel_top <= 0)
    {
        std::cerr << "print all features" << std::endl;
    }
    else if (sel_top < 1)
    {
        std::cerr << static_cast<int>(sel_top * 100 + 0.5) << "%" << std::endl;
    }
    else
    {
        std::cerr << static_cast<int>(sel_top + 0.5) << std::endl;
    }
    std::cerr << "Output:                       " << (out_path.empty() ? "to screen" : out_path) << ", ";
    std::cerr << (with_counts ? "with" : "without") << " count vectors" << std::endl
              << std::endl;
}

inline void ParseOptions(int argc, char *argv[],
                         std::string &idx_dir,
                         std::string &rk_mthd, size_t &nfold,
                         std::string &with_path, std::string &count_mode,
                         std::string &dsgn_path,
                         bool &ln_transf, bool &standardize, bool &no_norm,
                         float &sel_top,
                         std::string &out_path, bool &with_counts)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintRankHelper();
        exit(EXIT_SUCCESS);
    }
    size_t split_pos;
    std::string arg;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        arg = argv[i_opt];
        if (arg == "-h" || arg == "-help")
        {
            PrintRankHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-idxdir" && i_opt + 1 < argc)
        {
            idx_dir = argv[++i_opt];
        }
        else if (arg == "-rankby" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                nfold = std::stoi(arg.substr(split_pos + 1));
            }
            else
            {
                nfold = 1; // by default, without cross-validation
            }
            rk_mthd = arg.substr(0, split_pos);
        }
        else if (arg == "-with" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                count_mode = arg.substr(split_pos + 1);
            }
            else
            {
                count_mode = "rep"; // by default, take representative count vector (works also for general feature mode)
            }
            with_path = arg.substr(0, split_pos);
        }
        else if (arg == "-design" && i_opt + 1 < argc)
        {
            dsgn_path = argv[++i_opt];
        }
        else if (arg == "-ln")
        {
            ln_transf = true;
        }
        else if (arg == "-standardize")
        {
            standardize = true;
        }
        else if (arg == "-rankonraw")
        {
            no_norm = true;
        }
        else if (arg == "-seltop" && i_opt + 1 < argc)
        {
            sel_top = std::stof(argv[++i_opt]);
        }
        else if (arg == "-outpath" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else if (arg == "-withcounts")
        {
            with_counts = true;
        }
        else
        {
            PrintRankHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintRankHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (idx_dir.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    if (rk_mthd.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("-rankby STR is mandatory");
    }
    if (dsgn_path.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("-design STR is mandatory");
    }
}

#endif //KAMRAT_RANK_RANKRUNINFO_HPP
