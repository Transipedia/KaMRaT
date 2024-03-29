#ifndef KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP

void RankWelcome()
{
    std::cerr << "KaMRaT score: score features according to their association with sample conditions" << std::endl
              << "-----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintRankHelper()
{
    std::cerr << "[USAGE]    kamrat score -idxdir STR -count-mode STR -scoreby STR -design STR [-with STR1[:STR2] -seltop NUM -outpath STR -withcounts]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help             Print the helper" << std::endl;
    std::cerr << "            -idxdir STR          Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -scoreby STR         Scoring method, mandatory, can be one of: " << std::endl
              << "                                     classification (binary sample labels given by design file)" << std::endl
              << "                                         ttest.padj      adjusted p-value of t-test between conditions" << std::endl
              << "                                         ttest.pi        \u03C0-value of t-test between conditions" << std::endl
              << "                                         snr             signal-to-noise ratio between conditions" << std::endl
              << "                                         lr:nfold        accuracy by logistic regression classifier" << std::endl
              << "                                     classification (binary or multiple sample labels given by design file)" << std::endl
              << "                                         dids            DIDS score" << std::endl
              << "                                         bayes:nfold     accuracy by naive Bayes classifier" << std::endl
              << "                                     correlation evaluation (continuous sample labels given by design file)" << std::endl
              << "                                         pearson         Pearson correlation with the continunous sample condition" << std::endl
              << "                                         spearman        Spearman correlation with the continuous sample condition" << std::endl
              << "                                     unsupervised evaluation (no design file required)" << std::endl
              << "                                         sd              standard deviation" << std::endl
              << "                                         rsd1            standard deviation adjusted by mean" << std::endl
              << "                                         rsd2            standard deviation adjusted by min" << std::endl
              << "                                         rsd3            standard deviation adjusted by median" << std::endl
              << "                                         entropy         entropy of sample counts + 1" << std::endl;
    std::cerr << "            -design STR          Path to file indicating sample-condition design, mandatory unless using sd, rsd1, rsd2, rsd3, entropy" << std::endl
              << "                                     without header line, each row can be either: " << std::endl
              << "                                         sample name, sample condition" << std::endl
              << "                                         sample name, sample condition, sample batch (only for lrc, nbc, and svm)" << std::endl;
    std::cerr << "            -with STR1[:STR2]    File indicating features to score (STR1) and counting mode (STR2)" << std::endl
              << "                                     if not provided, all indexed features are used for scoring" << std::endl
              << "                                     STR2 can be one of [rep, mean, median]" << std::endl;
    std::cerr << "            -seltop NUM          Select top scored features" << std::endl
              << "                                     if NUM > 1, number of top features to select (should be integer)" << std::endl
              << "                                     if 0 < NUM <= 1, ratio of top features to select" << std::endl
              << "                                     if absent or NUM <= 0, output all features" << std::endl;
    std::cerr << "            -outpath STR         Path to scoring result" << std::endl
              << "                                     if not provided, output to screen" << std::endl;
    std::cerr << "            -withcounts          Output sample count vectors [false]" << std::endl
              << std::endl;
    std::cerr << "[NOTE]      For scoring methods lrc, nbc, and svm, a univariate CV fold number (nfold) can be provided" << std::endl
              << "                if nfold = 0, leave-one-out cross-validation" << std::endl
              << "                if nfold = 1, without cross-validation, training and testing on the whole datset" << std::endl
              << "                if nfold > 1, n-fold cross-validation" << std::endl
              << "            For t-test scoring methods, a transformation log2(x + 1) is applied to sample counts" << std::endl
              << "            For SVM scoring, sample counts standardization is applied feature by feature" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir,
                  const std::string &rk_mthd, const size_t nfold,
                  const std::string &with_path, const std::string &count_mode,
                  const std::string &dsgn_path,
                  const float sel_top,
                  const std::string &out_path, const bool with_counts)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                 " << idx_dir << std::endl;
    std::cerr << "Scoring method:               " << rk_mthd;
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

    std::cerr << "Scoring with:                 " << (with_path.empty() ? "features in index" : "features in " + with_path) << std::endl;
    std::cerr << "Feature counting mode:        " + count_mode << std::endl;
    if (!dsgn_path.empty())
    {
        std::cerr << "Sample design:                " << dsgn_path << std::endl;
    }
    std::cerr << "Selection of top features:    ";
    if (sel_top <= 0)
    {
        std::cerr << "print all features" << std::endl;
    }
    else if (sel_top < 1)
    {
        std::cerr << sel_top * 100 << "%" << std::endl;
    }
    else
    {
        std::cerr << static_cast<int>(sel_top + 0.5) << std::endl;
    }
    std::cerr << "Output:                       " << (out_path.empty() ? "to screen" : out_path) << ", ";
    std::cerr << (with_counts ? "with" : "without") << " count vectors" << std::endl
              << std::endl;
}

void ParseOptions(int argc, char *argv[],
                  std::string &idx_dir,
                  std::string &rk_mthd, size_t &nfold,
                  std::string &with_path, std::string &count_mode,
                  std::string &dsgn_path,
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
        else if (arg == "-scoreby" && i_opt + 1 < argc)
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
            with_path = arg.substr(0, split_pos);
        }
        else if (arg == "-design" && i_opt + 1 < argc)
        {
            dsgn_path = argv[++i_opt];
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
        throw std::invalid_argument("-scoreby STR is mandatory");
    }
    if (rk_mthd != "sd" && rk_mthd != "rsd1" && rk_mthd != "rsd2" && rk_mthd != "rsd3" && rk_mthd != "entropy" && dsgn_path.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("-design STR is mandatory");
    }
}

#endif //KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP
