#ifndef KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP

void RankWelcome()
{
    std::cerr << "KaMRaT rank: rank features according to their association with sample conditions" << std::endl
              << "-----------------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintRankHelper()
{
    std::cerr << "[USAGE]    kamrat rank -idxdir STR -count-mode STR -rankby STR -design STR [-with STR1[:STR2] -seltop NUM -outpath STR -withcounts]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help             Print the helper " << std::endl;
    std::cerr << "            -idxdir STR          Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -rankby STR          Ranking method, mandatory, can be one of: " << std::endl
              << "                                     ttest.padj      adjusted p-value of t-test between conditions" << std::endl
              << "                                     ttest.pi        \u03C0-value of t-test between conditions" << std::endl
              << "                                     snr             signal-to-noise ratio between conditions" << std::endl
              << "                                     dids            DIDS score" << std::endl
              << "                                     lr:nfold        accuracy by logistic regression classifier" << std::endl
              << "                                     bayes:nfold     accuracy by naive Bayes classifier" << std::endl
              << "                                     svm:nfold       accuracy on SVM classifier" << std::endl;
    std::cerr << "            -design STR          Path to file indicating sample-condition design" << std::endl
              << "                                     without header line, each row can be either: " << std::endl
              << "                                         sample name, sample condition" << std::endl
              << "                                         sample name, sample condition, sample batch (only for lrc, nbc, and svm)" << std::endl;
    std::cerr << "            -with STR1[:STR2]    File indicating features to rank (STR1) and counting mode (STR2)" << std::endl
              << "                                     if not provided, all indexed features are used for ranking" << std::endl
              << "                                     STR2 can be one of [rep, mean, median]" << std::endl;
    std::cerr << "            -seltop NUM          Select top ranked features" << std::endl
              << "                                     if NUM > 1, number of top features to select (should be integer)" << std::endl
              << "                                     if 0 < NUM <= 1, ratio of top features to select" << std::endl
              << "                                     if absent or NUM <= 0, output all features" << std::endl;
    std::cerr << "            -outpath STR         Path to ranking result" << std::endl
              << "                                     if not provided, output to screen" << std::endl;
    std::cerr << "            -withcounts          Output sample count vectors [false]" << std::endl
              << std::endl;
    std::cerr << "[NOTE]      For ranking methods lrc, nbc, and svm, a univariate CV fold number (nfold) can be provided" << std::endl
              << "                if nfold = 0, leave-one-out cross-validation" << std::endl
              << "                if nfold = 1, without cross-validation, training and testing on the whole datset" << std::endl
              << "                if nfold > 1, n-fold cross-validation" << std::endl
              << "            For t-test ranking methods, a transformation log2(x + 1) is applied to sample counts" << std::endl
              << "            For SVM ranking, sample counts standardization is applied feature by feature" << std::endl
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
        throw std::invalid_argument("-rankby STR is mandatory");
    }
    if (dsgn_path.empty())
    {
        PrintRankHelper();
        throw std::invalid_argument("-design STR is mandatory");
    }
}

#endif //KAMRAT_RUNINFOFILES_RANKRUNINFO_HPP
