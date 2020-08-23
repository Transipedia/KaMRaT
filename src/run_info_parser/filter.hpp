#ifndef KAMRAT_RUNINFOPARSER_RANK_HPP
#define KAMRAT_RUNINFOPARSER_RANK_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintFilterHelper()
{
    std::cerr << "=====> kamratFilter Helper <=====" << std::endl;
    std::cerr << "[USAGE]    kamratFilter -express-name STR -silent-name STR -express-thres INT_REC:INT_ABD -silent-thres INT_REC:INT_ABD [-smp-info STR] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]   -express-name all/rest/STR:cond/STR:smp   String indicating samples considered for express filter" << std::endl
              << "                                                         all         for considering all samples" << std::endl
              << "                                                         rest        for considering other samples not related with -silent-name" << std::endl
              << "                                                         STR:cond    for considering samples in the condition indicated by STR" << std::endl
              << "                                                         STR:smp     for considering the only sample indicated by STR" << std::endl;
    std::cerr << "           -silent-name all/rest/STR:cond/STR:smp    String indicating samples considered for silent filter" << std::endl
              << "                                                         all         for considering all samples" << std::endl
              << "                                                         rest        for considering other samples not related with -express-name" << std::endl
              << "                                                         STR:cond    for considering samples in the condition indicated by STR" << std::endl
              << "                                                         STR:smp     for considering the only one sample indicated by STR" << std::endl;
    std::cerr << "           -express-thres INT_REC:INT_ABD            Expressed filter threshold" << std::endl
              << "                                                         keep the feature with at least (>=) INT_REC samples labeled as express-name whose count >= INT_ABD" << std::endl
              << "                                                         if express-name indicates a sample name, INT_REC must be 1" << std::endl;
    std::cerr << "           -silent-thres INT_REC:INT_ABD             Silent filter threshold" << std::endl
              << "                                                         keep the feature with at least (>=) INT_REC samples labeled as silent-name whose count <= INT_ABD" << std::endl
              << "                                                         if silent-name indicates a sample name, INT_REC must be 1" << std::endl;
    std::cerr << "           -smp-info STR                             Path to sample-condition file, without header line" << std::endl
              << "                                                         if abscent, all columns except the first are regarded as samples and labeled as \"all\"" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &count_tab_path,
                         const std::string &express_level,
                         const std::string &express_name,
                         const int express_min_rec,
                         const int express_min_abd,
                         const std::string &silent_level,
                         const std::string &silent_name,
                         const int silent_min_rec,
                         const int silent_max_abd,
                         const std::string &sample_info_path)
{
    std::cerr << "k-mer count path:                           " << count_tab_path << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample info path:                           " << sample_info_path << std::endl;
    }
    std::cerr << "Express level:                              " << express_level << std::endl;
    std::cerr << "Express name:                               " << express_name << std::endl;
    std::cerr << "\tat least (>=) " << express_min_rec << " samples with count >= " << express_min_abd << std::endl;
    std::cerr << "Silent level:                               " << silent_level << std::endl;
    std::cerr << "Silent name:                                " << silent_name << std::endl;
    std::cerr << "\tat least (>=) " << silent_min_rec << " samples with count <= " << silent_max_abd << std::endl;
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &count_tab_path,
                         std::string &express_level,
                         std::string &express_name,
                         int &express_min_rec,
                         int &express_min_abd,
                         std::string &silent_level,
                         std::string &silent_name,
                         int &silent_min_rec,
                         int &silent_max_abd,
                         std::string &sample_info_path)
{
    int i_opt(1);
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintFilterHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-express-name" && i_opt + 1 < argc)
        {
            express_name = argv[++i_opt];
            if (express_name == "all" || express_name == "rest")
            {
                express_level = "global";
            }
            else
            {
                SubCommandParser(express_name, express_level);
                if (express_level != "cond" && express_level != "smp")
                {
                    PrintFilterHelper();
                    throw std::invalid_argument("express-name argument should have \'cond\' or \'smp\' followed by \':\' when filter by condition or sample");
                }
            }
        }
        else if (arg == "-silent-name" && i_opt + 1 < argc)
        {
            silent_name = argv[++i_opt];
            if (silent_name == "all" || silent_name == "rest")
            {
                silent_level = "global";
            }
            else
            {
                SubCommandParser(silent_name, silent_level);
                if (silent_level != "cond" && silent_level != "smp")
                {
                    PrintFilterHelper();
                    throw std::invalid_argument("silent-name argument should have \'cond\' or \'smp\' followed by \':\' when filter by condition or sample");
                }
            }
        }
        else if (arg == "-express-thres" && i_opt + 1 < argc)
        {
            std::string thres_rec = argv[++i_opt], thres_abd;
            SubCommandParser(thres_rec, thres_abd);
            express_min_rec = std::stoi(thres_rec);
            express_min_abd = std::stoi(thres_abd);
        }
        else if (arg == "-silent-thres" && i_opt + 1 < argc)
        {
            std::string thres_rec = argv[++i_opt], thres_abd;
            SubCommandParser(thres_rec, thres_abd);
            silent_min_rec = std::stoi(thres_rec);
            silent_max_abd = std::stoi(thres_abd);
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
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
    if (express_name.empty())
    {
        PrintFilterHelper();
        throw std::invalid_argument("express-name argument not specified");
    }
    if (express_min_rec == -1 || express_min_abd == -1)
    {
        PrintFilterHelper();
        throw std::invalid_argument("express-thres argument not specified or failed to parse");
    }
    if (silent_name.empty())
    {
        PrintFilterHelper();
        throw std::invalid_argument("silent-name argument not specified");
    }
    if (silent_min_rec == -1 || silent_max_abd == -1)
    {
        PrintFilterHelper();
        throw std::invalid_argument("silent-thres argument not specified or failed to parse");
    }
    if ((express_name == "all" && silent_name == "rest") || (silent_name == "all" && express_name == "rest"))
    {
        PrintFilterHelper();
        throw std::invalid_argument("all and rest are not compatible for express-/silent-names");
    }
    if ((express_level == "smp" && express_min_rec > 1) || (silent_level == "smp" && silent_min_rec > 1))
    {
        PrintFilterHelper();
        throw std::invalid_argument("Filtering in sample level but with min_recurrence > 1 does not make sense");
    }
}

#endif //KAMRAT_RUNINFOPARSER_RANK_HPP
