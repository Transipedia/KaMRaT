#ifndef KAMRAT_MERGE_MERGERUNINFO_HPP
#define KAMRAT_MERGE_MERGERUNINFO_HPP

#include <iostream>
#include <string>
#include <set>

#define MAX_PEARSON_DEFAULT 0.20
#define MAX_SPEARMAN_DEFAULT 0.26
#define MAX_MAC_DEFAULT 0.25

const std::set<std::string> kIntervMethodUniv{"pearson", "spearman", "mac"};
const std::set<std::string> kQuantModeUniv{"rep", "mean"};

const void PrintMergeHelper()
{
    std::cerr << "[Usage]    kamrat merge -klen INT -idx-path STR [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help               Print the helper" << std::endl;
    std::cerr << "            -klen INT              k-mer length (max_value: 32)" << std::endl;
    std::cerr << "            -idx-path STR          Temporary file path for saving count index, mandatory" << std::endl;
    std::cerr << "            -unstrand              If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -min-overlap INT       Min assembly overlap (max_value: k) [floor(k/2)]" << std::endl;
    std::cerr << "            -smp-info STR          Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                                       if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -interv-method STR     Intervention method (none, pearson, spearman, mac) [none]" << std::endl
              << "                                       the threshold can be precised after a ':' symbol" << std::endl;
    std::cerr << "            -quant-mode STR        Quantification mode (rep, mean) [rep]" << std::endl;
    std::cerr << "            -rep-name STR[:STR]    Representative value column name" << std::endl
              << "                                       if absent, k-mer input order will be taken as rep-val" << std::endl
              << "                                       if present, a representative mode can follow after a colon" << std::endl
              << "                                           min       take k-mer with minimum value as representative k-mer" << std::endl
              << "                                           minabs    take k-mer with minimum absolute value as representative k-mer" << std::endl
              << "                                           max       take k-mer with maximum value as representative k-mer" << std::endl
              << "                                           maxabs    take k-mer with maximum absolute value as representative k-mer" << std::endl;
    std::cerr << "            -out-path STR          Output contig count table path [default: output to screen]" << std::endl
              << std::endl;
}

const void PrintRunInfo(const size_t k_len,
                        const std::string &idx_path,
                        const bool stranded,
                        const size_t min_overlap,
                        const std::string &smp_info_path,
                        const std::string &interv_method, const float interv_thres,
                        const std::string &quant_mode,
                        const std::string &rep_colname,
                        const std::string &out_path,
                        const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:                  " << k_len << std::endl;
    std::cerr << "k-mer count index path:        " << idx_path << std::endl;
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

const void ParseOptions(int argc, char *argv[],
                        size_t &k_len,
                        std::string &idx_path,
                        bool &stranded,
                        size_t &min_overlap,
                        std::string &smp_info_path,
                        std::string &interv_method, float &interv_thres,
                        std::string &quant_mode,
                        std::string &rep_colname,
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
        }
        else if (arg == "-idx-path" && i_opt + 1 < argc)
        {
            idx_path = argv[++i_opt];
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
    if (idx_path.empty())
    {
        PrintMergeHelper();
        throw std::invalid_argument("temporary index path is mandatory");
    }

    if (min_overlap == 0)
    {
        min_overlap = (k_len / 2);
    }
    // dealing intervention method //
    std::string threshold_str;
    size_t split_pos = interv_method.find(":");
    if (split_pos != std::string::npos)
    {
        threshold_str = interv_method.substr(split_pos + 1);
        interv_method = interv_method.substr(0, split_pos);
    }
    if (interv_method != "none" && kIntervMethodUniv.find(interv_method) == kIntervMethodUniv.cend())
    {
        PrintMergeHelper();
        throw std::invalid_argument("unknown intervention method " + interv_method);
    }
    if (!rep_colname.empty() && smp_info_path.empty())
    {
        throw std::invalid_argument("when given representative column name, the sample-info path is mandatory");
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

#endif //KAMRAT_MERGE_MERGERUNINFO_HPP
