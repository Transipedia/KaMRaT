#ifndef KAMRAT_RUNINFOPARSER_MERGE_HPP
#define KAMRAT_RUNINFOPARSER_MERGE_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintMergeHelper()
{
    std::cerr << "[Usage]    kamrat merge -klen INT [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help              Print the helper" << std::endl;
    std::cerr << "            -klen INT             k-mer length (max_value: 32)" << std::endl;
    std::cerr << "            -unstrand             If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -min-overlap INT      Min assembly overlap (max_value: k) [floor(k/2)]" << std::endl;
    std::cerr << "            -smp-info STR         Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                                      if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -interv-method STR    Intervention method (none, pearson, spearman, mac) [none]" << std::endl
              << "                                      the threshold can be precised after a ':' symbol" << std::endl;
    std::cerr << "            -quant-mode STR       Quantification mode (rep, mean) [rep]" << std::endl;
    std::cerr << "            -rep-name STR         Representative value column name, k-mer input order as rep-val by default" << std::endl;
    std::cerr << "            -disk STR             Query on disk, followed by index file path" << std::endl;
    std::cerr << "            -out-path STR         Output contig count table path [default: output to screen]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const size_t k_len,
                         const bool stranded,
                         const size_t min_overlap,
                         const std::string &smp_info_path,
                         const std::string &interv_method, const float interv_thres,
                         const std::string &quant_mode,
                         const std::string &rep_colname,
                         const std::string &idx_path,
                         const std::string &out_path,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:                  " << k_len << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Min overlap length:            " << min_overlap << std::endl;
    if (!smp_info_path.empty())
    {
        std::cerr << "Sample-info path:              " << smp_info_path << std::endl;
    }
    std::cerr << "Intervention method:           " << interv_method << std::endl;
    if (interv_method != "none")
    {
        std::cerr << "\tintervention threshold = " << interv_thres << std::endl;
    }
    std::cerr << "Quantification mode:           " << quant_mode << std::endl;
    if (!rep_colname.empty())
    {
        std::cerr << "Representative k-mer:          with lowest value in column " << rep_colname << std::endl;
    }
    else
    {
        std::cerr << "Representative k-mer:          firstly inputted" << std::endl;
    }
    if (!idx_path.empty())
    {
        std::cerr << "Disk mode" << std::endl
                  << "\tindex path = " << idx_path << std::endl;
    }
    if (!out_path.empty())
    {
        std::cerr << "Output path:                   " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << "k-mer count path:              " << kmer_count_path << std::endl;
}

inline void ParseOptions(int argc, char *argv[],
                         size_t &k_len,
                         bool &stranded,
                         size_t &min_overlap,
                         std::string &smp_info_path,
                         std::string &interv_method, float &interv_thres,
                         std::string &quant_mode,
                         std::string &rep_colname,
                         std::string &idx_path,
                         std::string &out_path,
                         std::string &kmer_count_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h" || arg == "-help")
        {
            PrintMergeHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_len = atoi(argv[++i_opt]);
            if (min_overlap == 0)
            {
                min_overlap = (k_len / 2);
            }
        }
        else if (arg == "-unstrand")
        {
            stranded = false;
        }
        else if (arg == "-min-overlap" && i_opt + 1 < argc)
        {
            min_overlap = atoi(argv[++i_opt]);
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            smp_info_path = argv[++i_opt];
        }
        else if (arg == "-interv-method" && i_opt + 1 < argc)
        {
            interv_method = argv[++i_opt];
        }
        else if (arg == "-quant-mode" && i_opt + 1 < argc)
        {
            quant_mode = argv[++i_opt];
        }
        else if (arg == "-rep-name" && i_opt + 1 < argc)
        {
            rep_colname = argv[++i_opt];
        }
        else if (arg == "-disk" && i_opt + 1 < argc)
        {
            idx_path = argv[++i_opt];
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
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
    if (i_opt == argc)
    {
        PrintMergeHelper();
        throw std::invalid_argument("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];
    if (k_len == 0)
    {
        PrintMergeHelper();
        throw std::invalid_argument("k-mer length is mandatory");
    }

    // dealing intervention method //
    std::string threshold_str;
    SubCommandParser(interv_method, threshold_str);
    if (interv_method != "none" && INTERV_METHOD_UNIV.find(interv_method) == INTERV_METHOD_UNIV.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown intervention method " + interv_method);
    }
    if (threshold_str.empty() && interv_method == "pearson")
    {
        interv_thres = MAX_PEARSON_DEFAULT;
    }
    else if (threshold_str.empty() && interv_method == "spearman")
    {
        interv_thres = MAX_SPEARMAN_DEFAULT;
    }
    else if (threshold_str.empty() && interv_method == "mac")
    {
        interv_thres = MAX_MAC_DEFAULT;
    }
    else if (!threshold_str.empty())
    {
        interv_thres = std::stod(threshold_str);
    }
}

#endif //KAMRAT_RUNINFOPARSER_MERGE_HPP
