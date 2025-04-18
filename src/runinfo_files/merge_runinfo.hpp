#ifndef KAMRAT_RUNINFOFILES_MERGERUNINFO_HPP
#define KAMRAT_RUNINFOFILES_MERGERUNINFO_HPP

#include <unordered_set>

const std::unordered_set<std::string> kIntervMethodUniv{"none", "pearson", "spearman", "mac"};
const std::unordered_set<std::string> kRepModeUniv{"min", "minabs", "max", "maxabs"};
const std::unordered_set<std::string> kOutFmtUniv{"tab", "fa", "bin"};
const std::unordered_set<std::string> kCountModeUniv{"rep", "mean", "median"};
const std::unordered_set<std::string> kValueModeUniv{"int", "float"};

void MergeWelcome()
{
    std::cerr << "KaMRaT merge: extend k-mers into contigs" << std::endl
              << "------------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintMergeHelper()
{
    std::cerr << "[USAGE]    kamrat merge -idxdir STR -overlap MAX-MIN [-with STR1[:STR2] -interv STR[:FLOAT] -min-nbkmer INT -outfmt STR -outpath STR -counts STR1[:STR2]]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help               Print the helper" << std::endl;
    std::cerr << "            -idxdir STR            Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -overlap MAX-MIN       Overlap range for extension, by default: from (k-1) to \u230Ak/2\u230B" << std::endl
              << "                                       MIN and MAX are integers, MIN <= MAX < k-mer length" << std::endl;
    std::cerr << "            -with STR1[:STR2]      File indicating k-mers to be extended (STR1) and rep-mode (STR2)" << std::endl
              << "                                       if not provided, all indexed k-mers are used for extension" << std::endl
              << "                                       in the file STR1, a supplementary column of rep-value can be provided" << std::endl
              << "                                       STR2 can be one of {min, minabs, max, maxabs} [min]" << std::endl;
    std::cerr << "            -interv STR[:FLOAT]    Intervention method for extension [pearson:0.20]" << std::endl
              << "                                       can be one of {none, pearson, spearman, mac}" << std::endl
              << "                                       the threshold may follow a ':' symbol" << std::endl;
    std::cerr << "            -min-nbkmer INT        Minimal length of extended contigs [0]" << std::endl;
    std::cerr << "            -outfmt STR            Output format, STR can be `tab`, `fa`, or `bin` [default `tab`]" << std::endl
              << "                                       `tab` will output the final count table, set by default" << std::endl
              << "                                       `fa` will output a fasta file containing sequences without counts" << std::endl
              << "                                       `bin` will output a binary file, to be taken by the `-with` option of other modules" << std::endl;
    std::cerr << "            -outpath STR           Path to extension results" << std::endl
              << "                                       if not provided, output to screen" << std::endl;
    std::cerr << "            -counts STR1[:STR2]    How to compute contig counts from k-mer counts, only works if `-outfmt tab`" << std::endl
              << "                                       STR1 can be `rep`, `mean`, or `median` [default `rep`]" << std::endl
              << "                                           `rep` uses the representative k-mer count for the contig (see rep-mode in `-with`)" << std::endl
              << "                                           `mean` computes mean counts among all composite k-mers for each sample" << std::endl
              << "                                           `median` computes median counts among all composite k-mers for each sample" << std::endl
              << "                                       STR2 can be `int` or `float` [default `int`]" << std::endl
              << "                                           `int` will round the count values to nearest integers" << std::endl
              << "                                           `float` will output the values in decimals" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir, const size_t k_len, const bool stranded,
                  const size_t max_ovlp, const size_t min_ovlp,
                  const std::string &with_path, const std::string &rep_mode,
                  const std::string &itv_mthd, const float itv_thres,
                  const size_t min_nbkmer, const std::string &out_path,
                  const std::string &out_fmt, const std::string &count_mode, const std::string &value_mode)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                      " << idx_dir << std::endl;
    std::cerr << "k-mer length:                      " << k_len << std::endl;
    std::cerr << "Stranded extension:                " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Overlap range:                     from " << max_ovlp << " to " << min_ovlp << std::endl;
    std::cerr << "Merge k-mers in file:              " << (with_path.empty() ? "k-mers in index" : with_path) << std::endl;
    std::cerr << "Representative mode:               " << rep_mode << std::endl;
    std::cerr << "Intervention method:               " << itv_mthd
              << (itv_mthd != "none" ? (", threshold = " + std::to_string(itv_thres)) : "") << std::endl;
    std::cerr << "Minimal component k-mer number:    " + std::to_string(min_nbkmer) << std::endl;
    std::cerr << "Output:                            " << (out_path.empty() ? "to screen" : out_path) << std::endl
              << "    format:                        " + out_fmt << std::endl
              << "    count mode:                    " + count_mode << std::endl;
    if (out_fmt == "tab")
    {
        std::cerr << "    value mode:                    " + value_mode << std::endl;
    }
    std::cerr << std::endl;
}

void ParseOptions(int argc, char *argv[],
                  std::string &idx_dir, size_t &max_ovlp, size_t &min_ovlp,
                  std::string &with_path, std::string &rep_mode,
                  std::string &itv_mthd, float &itv_thres,
                  size_t &min_nbkmer, std::string &out_path,
                  std::string &out_fmt, std::string &count_mode, std::string &value_mode)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintMergeHelper();
        exit(EXIT_SUCCESS);
    }
    size_t split_pos;
    std::string arg;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        arg = argv[i_opt];
        if (arg == "-help" || arg == "-h")
        {
            PrintMergeHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-idxdir" && i_opt + 1 < argc)
        {
            idx_dir = argv[++i_opt];
        }
        else if (arg == "-overlap" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find("-");
            if (split_pos != std::string::npos)
            {
                min_ovlp = std::stoi(arg.substr(split_pos + 1));
                max_ovlp = std::stoi(arg.substr(0, split_pos));
            }
            else
            {
                throw std::invalid_argument("invalid overlap range: " + arg);
            }
            if (min_ovlp > max_ovlp)
            {
                throw std::invalid_argument("invalid overlap range, MAX should come first: " + arg);
            }
            if (min_ovlp > 31 || max_ovlp > 31)
            {
                throw std::invalid_argument("overlap range should not exceed 31");
            }
        }
        else if (arg == "-with" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                rep_mode = arg.substr(split_pos + 1);
            }
            with_path = arg.substr(0, split_pos);
        }
        else if (arg == "-interv" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                itv_thres = std::stof(arg.substr(split_pos + 1));
            }
            itv_mthd = arg.substr(0, split_pos);
        }
        else if (arg == "-min-nbkmer" && i_opt + 1 < argc)
        {
            min_nbkmer = std::stoul(argv[++i_opt]);
        }
        else if (arg == "-outfmt" && i_opt + 1 < argc)
        {
            out_fmt = argv[++i_opt];
        }
        else if (arg == "-outpath" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else if (arg == "-counts" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                value_mode = arg.substr(split_pos + 1);
            }
            count_mode = arg.substr(0, split_pos);
        }
        else
        {
            PrintMergeHelper();
            throw std::invalid_argument("unknown option " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintMergeHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (idx_dir.empty())
    {
        PrintMergeHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    // if (max_ovlp == 0 || min_ovlp == 0)
    // {
    //     PrintMergeHelper();
    //     throw std::invalid_argument("-overlap MAX-MIN is mandatory");
    // }
    if (kIntervMethodUniv.find(itv_mthd) == kIntervMethodUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown intervention method: " + itv_mthd);
    }
    if (kRepModeUniv.find(rep_mode) == kRepModeUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown representative mode: " + rep_mode);
    }
    if (kOutFmtUniv.find(out_fmt) == kOutFmtUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown output mode: " + out_fmt);
    }
    if (kCountModeUniv.find(count_mode) == kCountModeUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown count mode: " + count_mode);
    }
    if (kValueModeUniv.find(value_mode) == kValueModeUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown value mode: " + value_mode);
    }
}

#endif // KAMRAT_RUNINFOFILES_MERGERUNINFO_HPP
