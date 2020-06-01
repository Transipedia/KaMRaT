#ifndef KAMRAT_RUNINFOPARSER_EVALUATE_HPP
#define KAMRAT_RUNINFOPARSER_EVALUATE_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintHelper()
{
    std::cerr << "========= kamratEvaluate helper =========" << std::endl;
    std::cerr << "[Usage]    kamratEvaluate -l STRING -m STRING [-n] [-k INT] [-d STRING] STRING" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -l STRING  Contig list file path (MANDATORY)" << std::endl
              << "                       could be fasta, list, or table with contig sequences as the first column" << std::endl;
    std::cerr << "            -m STRING  Evaluate method:mode string (MANDATORY)" << std::endl
              << "                       method could be pearson, spearman, or mac" << std::endl
              << "                       mode could be farthest or worstAdj" << std::endl;
    std::cerr << "            -n         If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -k INT     k-mer length [31]" << std::endl;
    std::cerr << "            -d STRING  Colname list path indicating columns to print, either list or table with sample names as the first column" << std::endl
              << "                       if absent, all columns will be printed" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const std::string &contig_fasta_path,
                         const std::string &eval_method,
                         const std::string &eval_mode,
                         const bool stranded,
                         const unsigned int k_length,
                         const std::string &colname_list_path,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "Contig list path:       " << contig_fasta_path << std::endl;
    std::cerr << "Evaluation method:      " << eval_method << std::endl;
    std::cerr << "Evaluation mode:        " << eval_mode << std::endl;
    std::cerr << "Stranded mode:          " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:           " << k_length << std::endl;
    if (!colname_list_path.empty())
    {
        std::cerr << "Colname list path:      " << colname_list_path << std::endl;
    }
    std::cerr << "k-mer count path:       " << kmer_count_path << std::endl
              << std::endl;
}

inline bool ParseOptions(int argc,
                         char *argv[],
                         std::string &contig_fasta_path,
                         std::string &eval_method,
                         std::string &eval_mode,
                         bool &stranded,
                         unsigned int &k_length,
                         std::string &colname_list_path,
                         std::string &kmer_count_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h")
        {
            PrintHelper();
        }
        else if (arg == "-l" && i_opt + 1 < argc)
        {
            contig_fasta_path = argv[++i_opt];
        }
        else if (arg == "-m" && i_opt + 1 < argc)
        {
            eval_method = argv[++i_opt];
            SubCommandParser(eval_method, eval_mode);
        }
        else if (arg == "-n")
        {
            stranded = false;
        }
        else if (arg == "-k" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-d" && i_opt + 1 < argc)
        {
            colname_list_path = argv[++i_opt];
        }
        else
        {
            std::cerr << "ERROR: unknown option " << argv[i_opt] << std::endl;
            return false;
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        std::cerr << "ERROR: k-mer count table path is mandatory" << std::endl;
        return false;
    }
    kmer_count_path = argv[i_opt++];

    if (contig_fasta_path.empty())
    {
        std::cerr << "ERROR: contig path option is missing or failed to be parsed" << std::endl;
        return false;
    }
    if (EVALU_METHOD_UNIV.find(eval_method) == EVALU_METHOD_UNIV.cend())
    {
        std::cerr << "ERROR: evaluate method is missing or unknown" << std::endl;
        return false;
    }
    if (EVALU_MODE_UNIV.find(eval_mode) == EVALU_MODE_UNIV.cend())
    {
        std::cerr << "ERROR: evaluate mode is missing or unknown" << std::endl;
        return false;
    }
    return true;
}

#endif //KAMRAT_RUNINFOPARSER_EVALUATE_HPP
