#ifndef KAMRAT_MERGE_MERGERUNINFO_HPP
#define KAMRAT_MERGE_MERGERUNINFO_HPP

#include <iostream>
#include <string>
#include <unordered_set>

const std::unordered_set<std::string> kIntervMethodUniv{"none", "pearson", "spearman", "mac"};
const std::unordered_set<std::string> kSelectMode{"min", "minabs", "max", "maxabs"};

const void MergeWelcome()
{
    std::cerr << "KaMRaT merge: sequence extension" << std::endl
              << "----------------------------------------------------------------------------------------------------------" << std::endl;
}

const void PrintMergeHelper()
{
    std::cerr << "[Usage]    kamrat merge -overlap MAX-MIN -idx-dir STR [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help               Print the helper" << std::endl;
    std::cerr << "            -idxdir STR            Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -overlap MAX-MIN       Overlap range for extension, mandatory" << std::endl
              << "                                       MIN and MAX are integers, MIN <= MAX <= k-mer length" << std::endl;
    std::cerr << "            -selfile STR           File indicating k-mers to be extended" << std::endl
              << "                                       if absent, all indexed k-mers are used for extension" << std::endl
              << "                                       a second column with k-mer representativeness value can be provided" << std::endl;
    std::cerr << "            -repmode STR           Mode for selecting representative k-mers with value" << std::endl
              << "                                       can be one of: min, minabs, max, maxabs" << std::endl;
    std::cerr << "            -interv STR[:FLOAT]    Intervention method for extension [none]" << std::endl
              << "                                       could be one of (none, pearson, spearman, mac)" << std::endl
              << "                                       the threshold may follow a ':' symbol" << std::endl;
    std::cerr << "            -outpath STR           Path of extension result list" << std::endl
              << "                                       if absent, output to screen" << std::endl
              << std::endl;
}

const void PrintRunInfo(const std::string &idx_dir,
                        const size_t k_len,
                        const size_t max_ovlp, const size_t min_ovlp,
                        const bool stranded,
                        const std::string &sel_path, const std::string &rep_mode,
                        const std::string &itv_mthd, const float itv_thres,
                        const std::string &out_path)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:             " << idx_dir << std::endl;
    std::cerr << "k-mer length:             " << k_len << std::endl;
    std::cerr << "Overlap range:            from " << max_ovlp << " to " << min_ovlp << std::endl;
    std::cerr << "Stranded extension:       " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Select k-mers in file:    " << (sel_path.empty() ? "k-mers in index" : sel_path) << std::endl;
    std::cerr << "Intervention method:      " << itv_mthd << std::endl;
    if (itv_mthd != "none")
    {
        std::cerr << "\tthreshold = " << itv_thres << std::endl;
    }
    std::cerr << "Output:                 " << (out_path.empty() ? "to screen" : out_path) << std::endl;
}

const void ParseOptions(int argc, char *argv[],
                        std::string &idx_dir,
                        size_t &max_ovlp, size_t &min_ovlp,
                        std::string &sel_path, std::string &rep_mode,
                        std::string &itv_mthd, float &itv_thres,
                        std::string &out_path)
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
        else if (arg == "-selfile" && i_opt + 1 < argc)
        {
            sel_path = argv[++i_opt];
        }
        else if (arg == "-repmode" && i_opt + 1 < argc)
        {
            rep_mode = argv[++i_opt];
            if (kSelectMode.find(rep_mode) == kSelectMode.cend())
            {
                throw std::invalid_argument("select mode should be one of: min, minabs, max, maxabs");
            }
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
        else if (arg == "-outpath" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
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
        throw std::invalid_argument("KaMRaT index folder is mandatory");
    }
    if (max_ovlp == 0 || min_ovlp == 0)
    {
        PrintMergeHelper();
        throw std::invalid_argument("overlap range is mandatory");
    }
    if (kIntervMethodUniv.find(itv_mthd) == kIntervMethodUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown intervention method: " + itv_mthd);
    }
}

#endif //KAMRAT_MERGE_MERGERUNINFO_HPP
