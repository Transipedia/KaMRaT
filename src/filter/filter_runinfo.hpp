#ifndef KAMRAT_FILTER_FILTERRUNINFO_HPP
#define KAMRAT_FILTER_FILTERRUNINFO_HPP

void FilterWelcome()
{
    std::cerr << "KaMRaT filter: feature filter by expression level" << std::endl
              << "-----------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintFilterHelper()
{
    std::cerr << "[USAGE]    kamrat filter -filter-info STR [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help              Print the helper" << std::endl;
    std::cerr << "            -idxdir STR           Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -design STR           Path to filter design file, a table of two columns, mandatory" << std::endl
              << "                                      the first column indicate sample names" << std::endl
              << "                                      the second column should be either UP or DOWN (capital letters)" << std::endl
              << "                                          samples with UP will be considered as up-regulated samples" << std::endl
              << "                                          samples with DOWN will be considered as down-regulated samples" << std::endl
              << "                                          samples not given will be neutral (not considered for filter)" << std::endl
              << "                                          samples can also be all UP or all DOWN" << std::endl;
    std::cerr << "            -upmin INT1:INT2      Up feature lower bound, [1:1, meaning no filter]" << std::endl
              << "                                      output features counting >= INT1 in >= INT2 UP-samples" << std::endl;
    std::cerr << "            -downmax INT1:INT2    Down feature upper bound [inf:1, meaning no filter]" << std::endl
              << "                                      output features counting <= INT1 in >= INT2 DOWN-samples" << std::endl;
    std::cerr << "            -normalize            Filter with normalized counts [false]" << std::endl;
    std::cerr << "            -outpath STR          Path to results after filter" << std::endl
              << "                                      if not provided, output to screen" << std::endl;
    std::cerr << "            -withcounts           Output sample count vectors [false]" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir,
                        const std::string &dsgn_path,
                        const size_t up_min_abd, const size_t up_min_rec,
                        const size_t down_max_abd, const size_t down_min_rec,
                        const bool norm,
                        const std::string &out_path, const bool with_counts)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                        " << idx_dir << std::endl;
    std::cerr << "Path to filter design file:          " << dsgn_path << std::endl;
    std::cerr << "Up-regulated lower bound:            " << std::endl
              << "\tfeatures counting >= " << up_min_abd << " in >= " << up_min_rec << " up-regulated samples" << std::endl;
    std::cerr << "Down-regulated upper bound:          " << std::endl
              << "\tfeatures counting <= " << (down_max_abd == std::numeric_limits<size_t>::max() ? "inf" : std::to_string(down_max_abd))
              << " in >= " << down_min_rec << "down-regulated samples" << std::endl;
    std::cerr << "Filter after count normalization:    " << (norm ? "TRUE" : "FALSE") << std::endl;
    std::cerr << "Output:                              " << (out_path.empty() ? "to screen" : out_path) << ", ";
    std::cerr << (with_counts ? "with" : "without") << " count vectors" << std::endl
              << std::endl;
}

void ParseOptions(int argc, char *argv[],
                        std::string &idx_dir,
                        std::string &dsgn_path,
                        size_t &up_min_abd, size_t &up_min_rec,
                        size_t &down_max_abd, size_t &down_min_rec,
                        bool &norm,
                        std::string &out_path, bool &with_counts)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintFilterHelper();
        exit(EXIT_SUCCESS);
    }
    size_t split_pos;
    std::string arg;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        arg = argv[i_opt];
        if (arg == "-help" || arg == "-h")
        {
            PrintFilterHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-design" && i_opt + 1 < argc)
        {
            dsgn_path = argv[++i_opt];
        }
        else if (arg == "-upmin" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                up_min_abd = std::stoul(arg.substr(0, split_pos));
                up_min_rec = std::stoul(arg.substr(split_pos + 1));
            }
            else
            {
                throw std::invalid_argument("unable to parse -upmin argument: " + arg);
            }
        }
        else if (arg == "-downmax" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                down_max_abd = std::stoul(arg.substr(0, split_pos));
                down_min_rec = std::stoul(arg.substr(split_pos + 1));
            }
            else
            {
                throw std::invalid_argument("unable to parse -downmax argument: " + arg);
            }
        }
        else if (arg == "-norm")
        {
            norm = true;
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
            PrintFilterHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintFilterHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (idx_dir.empty())
    {
        PrintFilterHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    if (dsgn_path.empty())
    {
        PrintFilterHelper();
        throw std::domain_error("-design STR is mandatory");
    }
}

#endif //KAMRAT_FILTER_FILTERRUNINFO_HPP
