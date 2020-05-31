#ifndef KAMRAT_RUNINFOPARSER_NORM_HPP
#define KAMRAT_RUNINFOPARSER_NORM_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratNorm helper =========" << std::endl;
    std::cerr << "[Usage]    kamratNorm -b CHAR [-d STRING] [-T STRING] STRING" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -b CHAR    BASE for normalization: count per BASE (MANDATORY, could be B, M, or K)" << std::endl;
    std::cerr << "            -d STRING  Sample-info path" << std::endl
              << "                       if absent, all columns will be printed" << std::endl;
    std::cerr << "            -T STRING  Transformation before evaluation (e.g. log)" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const size_t &baseN,
                         const std::string &sample_info_path,
                         const std::string &trans_mode,
                         const std::string &count_tab_path)
{
    std::cerr << std::endl;
    std::cerr << "Normalization base:     " << baseN << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample-info path:      " << sample_info_path << std::endl;
    }
    if (!trans_mode.empty())
    {
        std::cerr << "Transformation mode:    " << trans_mode << std::endl;
    }
    std::cerr << "k-mer count path:       " << count_tab_path << std::endl
              << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &baseN,
                         std::string &sample_info_path,
                         std::string &trans_mode,
                         std::string &count_tab_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h")
        {
            PrintHelper();
        }
        else if (arg == "-b" && i_opt + 1 < argc)
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
        else if (arg == "-d" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
        }
        else if (arg == "-T" && i_opt + 1 < argc)
        {
            trans_mode = argv[++i_opt];
        }
        else
        {
            throw std::domain_error("unknown option " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];

    if (baseN == 0)
    {
        throw std::domain_error("base mode is missing or failed to be parsed");
    }
    if (!trans_mode.empty() && trans_mode != "log")
    {
        throw std::domain_error("unknown transformation mode " + trans_mode);
    }
}

#endif //KAMRAT_RUNINFOPARSER_NORM_HPP
