#ifndef KAMRAT_RUNINFOPARSER_NORM_HPP
#define KAMRAT_RUNINFOPARSER_NORM_HPP

#include <iostream>
#include <string>

inline void PrintNormHelper()
{
    std::cerr << "========= kamratNorm helper =========" << std::endl;
    std::cerr << "[Usage]    kamratNorm -base CHAR [-smp-info STR] [-ln] [-smp-sum STR] COUNT_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help       Print the helper" << std::endl;
    std::cerr << "            -base CHAR     BASE for normalization (MANDATORY)" << std::endl
              << "                               B: count per billion" << std::endl
              << "                               M: count per million" << std::endl
              << "                               K: count per thousand" << std::endl;
    std::cerr << "            -smp-info STR  Sample-info path" << std::endl
              << "                               if absent, all columns except the first will be normalized" << std::endl;
    std::cerr << "            -ln            Apply ln(x + 1) transformation after normalization [false]" << std::endl
              << "                               hint: please remember to unlog for original counts after kamratReduce or other analysis" << std::endl;
    std::cerr << "            -smp-sum STR   Path for outputing sample sum [./sample_sum.tsv]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const size_t &baseN,
                         const std::string &sample_info_path,
                         const bool ln_trans,
                         const std::string &sample_sum_path,
                         const std::string &count_tab_path)
{
    std::cerr << std::endl;
    std::cerr << "Normalization base:     " << baseN << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample-info path:       " << sample_info_path << std::endl;
    }
    std::cerr << "Ln transformation:      " << (ln_trans ? "On" : "Off") << std::endl;
    std::cerr << "Sample sum path:        " << sample_sum_path << std::endl;
    std::cerr << "k-mer count path:       " << count_tab_path << std::endl
              << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &baseN,
                         std::string &sample_info_path,
                         bool &ln_trans,
                         std::string &sample_sum_path,
                         std::string &count_tab_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h" || arg == "-help")
        {
            PrintNormHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-base" && i_opt + 1 < argc)
        {
            char base_mode = *(argv[++i_opt]);
            switch (base_mode)
            {
            case 'B':
                baseN = 1E9;
                break;
            case 'M':
                baseN = 1E6;
                break;
            case 'K':
                baseN = 1E3;
                break;
            default:
                throw std::domain_error("unknown base mode character " + base_mode);
            }
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
        }
        else if (arg == "-ln")
        {
            ln_trans = true;
        }
        else if (arg == "-smp-sum" && i_opt + 1 < argc)
        {
            sample_sum_path = argv[++i_opt];
        }
        else
        {
            PrintNormHelper();
            throw std::domain_error("unknown option " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintNormHelper();
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];

    if (baseN == 0)
    {
        PrintNormHelper();
        throw std::domain_error("base mode is missing or failed to be parsed");
    }
}

#endif //KAMRAT_RUNINFOPARSER_NORM_HPP
