#ifndef KAMRAT_FILTER_FILTERRUNINFO_HPP
#define KAMRAT_FILTER_FILTERRUNINFO_HPP

#include <iostream>
#include <string>

const void PrintFilterHelper()
{
    std::cerr << "[USAGE]    kamrat filter -filter-info STR [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]   -h,-help               Print the helper" << std::endl;
    std::cerr << "           -filter-info STR       Filter-info path, should be a table of two columns WITHOUT header row, MANDATORY" << std::endl
              << "                                      the first column should be sample names" << std::endl
              << "                                      the second column should be either UP or DOWN (all capital letters)" << std::endl
              << "                                          samples with UP will be considered as up-regulated samples" << std::endl
              << "                                          samples with DOWN will be considered as down-regulated samples" << std::endl
              << "                                          samples not mentioned in the file will not taken into consideration (being neutral)" << std::endl
              << "                                          samples can also be all UP or all DOWN" << std::endl;
    std::cerr << "           -up-min INT1:INT2      Up feature lower bound, [by default 1:1 (no filter)]" << std::endl
              << "                                      print the feature if at least (>=) INT1 UP-samples have count >= INT2" << std::endl;
    std::cerr << "           -down-max INT1:INT2    Down feature upper bound [by default 1:inf (no filter)]" << std::endl
              << "                                      print the feature if at least (>=) INT1 DOWN-samples have count <= INT2" << std::endl;
    std::cerr << "           -out-path              Output table path [default: output to screen]" << std::endl
              << std::endl;
}

const void PrintRunInfo(const std::string &filter_info_path,
                        const size_t up_min_rec, const size_t up_min_abd,
                        const size_t down_min_rec, const size_t down_max_abd,
                        const std::string &out_path,
                        const std::string &count_tab_path)
{
    std::cerr << "Filter-info path:                           " << filter_info_path << std::endl;
    std::cerr << "Up-regulated lower bound:                   " << std::endl
              << "    at least " << up_min_rec << " up-regulated samples are counted at least " << up_min_abd << std::endl;
    std::cerr << "Down-regulated upper bound:                 " << std::endl
              << "    at least " << down_min_rec << " down-regulated samples are counted at most " << down_max_abd << std::endl;
    if (!out_path.empty())
    {
        std::cerr << "Output path:                   " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << "Count table path:                           " << count_tab_path << std::endl;
    std::cerr << std::endl;
}

const void ParseOptions(int argc, char *argv[],
                        std::string &filter_info_path,
                        size_t &up_min_rec, size_t &up_min_abd,
                        size_t &down_min_rec, size_t &down_max_abd,
                        std::string &out_path,
                        std::string &count_tab_path)
{
    std::string str_thres;
    int i_opt(1);
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintFilterHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-filter-info" && i_opt + 1 < argc)
        {
            filter_info_path = argv[++i_opt];
        }
        else if (arg == "-up-min" && i_opt + 1 < argc)
        {
            str_thres = argv[++i_opt];
            size_t split_pos = str_thres.find(":");
            if (split_pos != std::string::npos)
            {
                up_min_rec = std::stoul(str_thres.substr(0, split_pos));
                up_min_abd = std::stoul(str_thres.substr(split_pos + 1));
            }
        }
        else if (arg == "-down-max" && i_opt + 1 < argc)
        {
            str_thres = argv[++i_opt];
            size_t split_pos = str_thres.find(":");
            if (split_pos != std::string::npos)
            {
                down_min_rec = std::stoul(str_thres.substr(0, split_pos));
                down_max_abd = std::stoul(str_thres.substr(split_pos + 1));
            }
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else
        {
            PrintFilterHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintFilterHelper();
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];
}

#endif //KAMRAT_FILTER_FILTERRUNINFO_HPP
