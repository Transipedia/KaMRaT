#ifndef KAMRAT_RUNINFOPARSER_MERGE_HPP
#define KAMRAT_RUNINFOPARSER_MERGE_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintMergeHelper()
{
    std::cerr << "[Usage]    kamrat merge -klen INT [-min-overlap INT] [-unstrand] [-smp-info STR] [-interv STR] [-quant STR] [-rep-name STR] [-disk] [-idx-dir STR] KMER_COUNT_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help            Print the helper" << std::endl;
    std::cerr << "            -klen INT           k-mer length (max_value: 32)" << std::endl;
    std::cerr << "            -unstrand           If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -min-overlap INT    Min assembly overlap (max_value: k) [floor(k/2)]" << std::endl;
    std::cerr << "            -smp-info STR       Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                                    if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -interv STR         Intervention method (none, pearson, spearman, mac) [none]" << std::endl
              << "                                    the threshold can be precised after a ':' symbol" << std::endl;
    std::cerr << "            -quant STR          Quantification mode (rep, mean) [rep]" << std::endl;
    std::cerr << "            -rep-name STR       Representative value column name, k-mer input order as rep-val by default" << std::endl;
    std::cerr << "            -disk               Query on disk [false]" << std::endl;
    std::cerr << "            -idx-dir STR        Count index directory path [./]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const size_t k_length,
                         const bool stranded,
                         const size_t min_overlap,
                         const std::string &sample_info_path,
                         const std::string &interv_method,
                         const float interv_thres,
                         const std::string &quant_mode,
                         const std::string &rep_value_cname,
                         const bool disk_mode,
                         const std::string &tmp_dir,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Min overlap length:            " << min_overlap << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample-info path:              " << sample_info_path << std::endl;
    }
    std::cerr << "Intervention method:           " << interv_method << std::endl;
    if (interv_method != "none")
    {
        std::cerr << "Intervention threshold:        " << interv_thres << std::endl;
    }
    std::cerr << "Quantification mode:           " << quant_mode << std::endl;
    if (!rep_value_cname.empty())
    {
        std::cerr << "Representative value:          " << rep_value_cname << std::endl;
    }
    else if (quant_mode != "mean")
    {
        std::cerr << "Representative value:          input order" << std::endl;
    }
    std::cerr << "Disk mode:                     " << (disk_mode ? "On" : "Off") << std::endl;
    if (!tmp_dir.empty())
    {
        std::cerr << "temporary directory path:      " << tmp_dir << std::endl;
    }
    std::cerr << "k-mer count path:              " << kmer_count_path << std::endl
              << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &k_length,
                         bool &stranded,
                         size_t &min_overlap,
                         std::string &sample_info_path,
                         std::string &interv_method,
                         float &interv_thres,
                         std::string &quant_mode,
                         std::string &rep_value_cname,
                         bool &disk_mode,
                         std::string &tmp_dir,
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
            k_length = atoi(argv[++i_opt]);
            if (min_overlap == 0)
            {
                min_overlap = (k_length / 2);
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
            sample_info_path = argv[++i_opt];
        }
        else if (arg == "-interv" && i_opt + 1 < argc)
        {
            interv_method = argv[++i_opt];
        }
        else if (arg == "-quant" && i_opt + 1 < argc)
        {
            quant_mode = argv[++i_opt];
        }
        else if (arg == "-rep-name" && i_opt + 1 < argc)
        {
            rep_value_cname = argv[++i_opt];
        }
        else if (arg == "-disk")
        {
            disk_mode = true;
        }
        else if (arg == "-idx-dir" && i_opt + 1 < argc)
        {
            tmp_dir = argv[++i_opt];
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
    if (k_length == 0)
    {
        PrintMergeHelper();
        throw std::invalid_argument("k-mer length is missing");
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
        interv_thres = MIN_PEARSON_DEFAULT;
    }
    else if (threshold_str.empty() && interv_method == "spearman")
    {
        interv_thres = MIN_SPEARMAN_DEFAULT;
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
