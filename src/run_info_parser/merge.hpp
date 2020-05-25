#ifndef KAMRAT_RUNINFOPARSER_MERGE_HPP
#define KAMRAT_RUNINFOPARSER_MERGE_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintHelper(const size_t k_length, const size_t min_overlap)
{
    std::cerr << "========= kamratMerge helper =========" << std::endl;
    std::cerr << "[Usage]    kamratMerge [-k k_length] [-m min_overlap] [-n] [-x] [-d sample_info] [-i interv_method] [-q quant_mode] [-t tmp_dir] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -n         If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -k INT     k-mer length (max_value: 32) [" << int(k_length) << "]" << std::endl;
    std::cerr << "            -m INT     Min assembly overlap (max_value: k) [" << int(min_overlap) << "]" << std::endl;
    std::cerr << "            -d STRING  Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                       if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -i STRING  Intervention method (none, pearson, spearman, mac) [none]" << std::endl
              << "                       the threshold can be precised after a ':' symbol" << std::endl;
    std::cerr << "            -q STRING  Quantification mode (rep, mean) [rep]" << std::endl;
    std::cerr << "            -x         Query on disk [false]" << std::endl;
    std::cerr << "            -t STRING  Temporary directory [./]" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const bool stranded,
                         const size_t k_length,
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
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
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
        std::cerr << "Representative column name:    " << rep_value_cname << std::endl;
    }
    else if (quant_mode != "mean")
    {
        std::cerr << "Representative column name:    input order" << std::endl;
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
                         bool &stranded,
                         size_t &k_length,
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
        if (arg == "-h")
        {
            PrintHelper(k_length, min_overlap);
        }
        else if (arg == "-n")
        {
            stranded = false;
        }
        else if (arg == "-k" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-m" && i_opt + 1 < argc)
        {
            min_overlap = atoi(argv[++i_opt]);
        }
        else if (arg == "-d" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
        }
        else if (arg == "-i" && i_opt + 1 < argc)
        {
            interv_method = argv[++i_opt];
        }
        else if (arg == "-q" && i_opt + 1 < argc)
        {
            quant_mode = argv[++i_opt];
        }
        else if (arg == "-x")
        {
            disk_mode = true;
        }
        else if (arg == "-t" && i_opt + 1 < argc)
        {
            tmp_dir = argv[++i_opt];
        }
        else
        {
            throw std::domain_error("unknown option " + arg);
            exit(EXIT_FAILURE);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        throw std::domain_error("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];

    // dealing quantification mode //
    SubCommandParser(quant_mode, rep_value_cname);
    if (QUANT_MODE_UNIV.find(quant_mode) == QUANT_MODE_UNIV.cend())
    {
        throw std::domain_error("unknown quant mode " + quant_mode);
    }

    // dealing intervention method //
    std::string threshold_str;
    SubCommandParser(interv_method, threshold_str);
    if (INTERV_METHOD_UNIV.find(interv_method) == INTERV_METHOD_UNIV.cend())
    {
        throw std::domain_error("unknown intervention method " + interv_method);
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