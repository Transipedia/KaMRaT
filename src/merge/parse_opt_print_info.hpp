#ifndef KAMRAT_MERGE_PARSEOPTPRINTINFO_HPP
#define KAMRAT_MERGE_PARSEOPTPRINTINFO_HPP

#include <iostream>
#include <string>

#include "utils.hpp"

inline void PrintHelper(const size_t k_length, const size_t min_overlap)
{
    std::cerr << "========= kamratMerge helper =========" << std::endl;
    std::cerr << "[Usage]    kamratMerge [-k k_length] [-m min_overlap] [-nj] [-d sample_info] [-i interv_method] [-q quant_mode] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -n         If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -k INT     k-mer length (max_value: 32) [" << int(k_length) << "]" << std::endl;
    std::cerr << "            -m INT     Min assembly overlap (max_value: k) [" << int(min_overlap) << "]" << std::endl;
    std::cerr << "            -d STRING  Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                       if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -i STRING  Intervention method (" << INTERV_NONE << ", " << INTERV_PEARSON << ", " << INTERV_SPEARMAN
              << ", " << INTERV_MAC << ") [" << INTERV_NONE << "]" << std::endl
              << "                       the threshold can be precised after a ':' symbol" << std::endl;
    std::cerr << "            -q STRING  Quantification mode (rep, mean) [rep]" << std::endl;
    std::cerr << "            -j         Adjacent k-mer comparison (valid only with intervention) [false]" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void SubCommandParser(std::string &command, std::string &sub_command)
{
    size_t split_pos = command.find(":");
    if (split_pos != std::string::npos)
    {
        sub_command = command.substr(split_pos + 1);
        command = command.substr(0, split_pos);
    }
}

inline void PrintRunInfo(const bool stranded,
                         const bool disk_mode,
                         const size_t k_length,
                         const size_t min_overlap,
                         const std::string &sample_info_path,
                         const std::string &tmp_dir,
                         const std::string &interv_method,
                         const double interv_thres,
                         const std::string &quant_mode,
                         const bool comp_adj,
                         const std::string &rep_value_cname,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Disk mode:                     " << (disk_mode ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
    std::cerr << "Min overlap length:            " << min_overlap << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample-info path:              " << sample_info_path << std::endl;
    }
    std::cerr << "Intervention method:           " << interv_method << std::endl;
    if (interv_method != INTERV_NONE)
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
    std::cerr << "Adjacent comparison:           " << (comp_adj ? "On" : "Off") << std::endl;
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
                         bool &disk_mode,
                         size_t &k_length,
                         size_t &min_overlap,
                         std::string &sample_info_path,
                         std::string &tmp_dir,
                         std::string &interv_method,
                         double &interv_thres,
                         std::string &quant_mode,
                         bool &comp_adj,
                         std::string &rep_value_cname,
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
        // else if (arg == "-x")
        // {
        //     disk_mode = true;
        // }
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
        else if (arg == "-j")
        {
            comp_adj = true;
        }
        // else if (arg == "-t" && i_opt + 1 < argc)
        // {
        //     tmp_dir = argv[++i_opt];
        // }
        else
        {
            std::cerr << "ERROR: unknown option " << argv[i_opt] << std::endl;
            exit(EXIT_FAILURE);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        std::cerr << "ERROR: k-mer count table path is mandatory" << std::endl;
        exit(EXIT_FAILURE);
    }
    kmer_count_path = argv[i_opt++];

    // dealing quantification mode //
    SubCommandParser(quant_mode, rep_value_cname);
    if (quant_mode != "rep" && quant_mode != "mean")
    {
        std::cerr << "ERROR: unknown quant mode " << quant_mode << std::endl;
        exit(EXIT_FAILURE);
    }

    // dealing intervention method //
    std::string threshold_str;
    SubCommandParser(interv_method, threshold_str);
    if (interv_method != INTERV_NONE && interv_method != INTERV_PEARSON && interv_method != INTERV_SPEARMAN && interv_method != INTERV_MAC)
    {
        std::cerr << "ERROR: unknown intervention method " << interv_method << std::endl;
        exit(EXIT_FAILURE);
    }
    if (threshold_str.empty() && interv_method == INTERV_PEARSON)
    {
        interv_thres = MIN_PEARSON;
    }
    else if (threshold_str.empty() && interv_method == INTERV_SPEARMAN)
    {
        interv_thres = MIN_SPEARMAN;
    }
    else if (threshold_str.empty() && interv_method == INTERV_MAC)
    {
        interv_thres = MAX_MAC;
    }
    else if (!threshold_str.empty())
    {
        interv_thres = std::stod(threshold_str);
    }
}

#endif //KAMRAT_MERGE_PARSEOPTPRINTINFO_HPP