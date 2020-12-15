#ifndef KAMRAT_RUNINFOPARSER_EVALUATE_HPP
#define KAMRAT_RUNINFOPARSER_EVALUATE_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintHelper()
{
    std::cerr << "========= kamratEvaluate helper =========" << std::endl;
    std::cerr << "[Usage]    contigEvaluate -fasta STR -eval-method STR -idx-path STR [-options] KMER_COUNT_TAB" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h, -help           Print the helper" << std::endl;
    std::cerr << "            -fasta STR          Contig fasta file path (MANDATORY)" << std::endl;
    std::cerr << "            -eval-method STR    Evaluate method:mode string (MANDATORY)" << std::endl
              << "                                    method could be pearson, spearman, or mac" << std::endl
              << "                                    mode could be farthest or worstAdj" << std::endl;
    std::cerr << "            -idx-path STR       Count index file path" << std::endl;
    std::cerr << "            -unstrand           If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -klen INT               k-mer length [31]" << std::endl;
    std::cerr << "            -smp-info STR       Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                                    if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -out-path STR       Output table path [default: output to screen]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &contig_fasta_path,
                         const std::string &eval_method,
                         const std::string &eval_mode,
                         const std::string &idx_path,
                         const bool stranded,
                         const unsigned int k_length,
                         const std::string &colname_list_path,
                         const std::string &out_path,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "Contig list path:        " << contig_fasta_path << std::endl;
    std::cerr << "Evaluation method:       " << eval_method << std::endl;
    std::cerr << "Evaluation mode:         " << eval_mode << std::endl;
    std::cerr << "Index directory path:    " << idx_path << std::endl;
    std::cerr << "Stranded mode:           " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:            " << k_length << std::endl;
    if (!colname_list_path.empty())
    {
        std::cerr << "Colname list path:       " << colname_list_path << std::endl;
    }
    if (!out_path.empty())
    {
        std::cerr << "Output path:             " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << "k-mer count path:        " << kmer_count_path << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &contig_fasta_path,
                         std::string &eval_method,
                         std::string &eval_mode,
                         std::string &idx_path,
                         bool &stranded,
                         unsigned int &k_length,
                         std::string &colname_list_path,
                         std::string &out_path,
                         std::string &kmer_count_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            contig_fasta_path = argv[++i_opt];
        }
        else if (arg == "-eval-method" && i_opt + 1 < argc)
        {
            eval_method = argv[++i_opt];
            SubCommandParser(eval_method, eval_mode);
        }
        else if (arg == "-idx-path" && i_opt + 1 < argc)
        {
            idx_path = argv[++i_opt];
        }
        else if (arg == "-unstrand")
        {
            stranded = false;
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            colname_list_path = argv[++i_opt];
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else
        {
            PrintHelper();
            throw std::invalid_argument("unknown option " + std::string(argv[i_opt]));
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintHelper();
        throw std::invalid_argument("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];

    if (contig_fasta_path.empty())
    {
        PrintHelper();
        throw std::invalid_argument("contig path option is missing or failed to be parsed");
    }
    if (EVALU_METHOD_UNIV.find(eval_method) == EVALU_METHOD_UNIV.cend())
    {
        PrintHelper();
        throw std::invalid_argument("evaluate method is missing or unknown");
    }
    if (EVALU_MODE_UNIV.find(eval_mode) == EVALU_MODE_UNIV.cend())
    {
        PrintHelper();
        throw std::invalid_argument("evaluate mode is missing or unknown");
    }
}

#endif //KAMRAT_RUNINFOPARSER_EVALUATE_HPP
