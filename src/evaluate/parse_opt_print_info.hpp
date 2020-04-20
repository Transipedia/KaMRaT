#ifndef KAMRAT_EVALUATE_PARSEOPTPRINTINFO_HPP
#define KAMRAT_EVALUATE_PARSEOPTPRINTINFO_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratEvaluate helper =========" << std::endl;
    std::cerr << "[Usage]    kamratEvaluate -l contig_list_path -m eval_method [-n] [-k k_length] [-d colname_list_path] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -l STRING  Contig list file path (MANDATORY)" << std::endl
              << "                       could be fasta, list, or table with contig sequences as the first column" << std::endl;
    std::cerr << "            -m STRING  Evaluate method (mean, adj-mac, adj-pearson, adj-spearman)" << std::endl;
    std::cerr << "            -n         If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -k INT     k-mer length [31]" << std::endl;
    std::cerr << "            -d STRING  Colname list path indicating columns to print, either list or table with sample names as the first column" << std::endl
              << "                       if absent, all columns will be printed" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const std::string &contig_list_path,
                         const unsigned int k_length,
                         const unsigned int min_overlap,
                         const bool stranded,
                         const std::string &colname_list_path,
                         const std::string &eval_method,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "Contig list path:       " << contig_list_path << std::endl;
    std::cerr << "Evaluate method:        " << eval_method << std::endl;
    std::cerr << "Stranded mode:          " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:           " << k_length << std::endl;
    std::cerr << "Minimum overlap:        " << min_overlap << std::endl;
    if (!colname_list_path.empty())
    {
        std::cerr << "Colname list path:      " << colname_list_path << std::endl;
    }
    std::cerr << "k-mer count path:       " << kmer_count_path << std::endl
              << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &contig_list_path,
                         std::string &eval_method,
                         bool &stranded,
                         unsigned int &k_length,
                         unsigned int &min_overlap,
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
            contig_list_path = argv[++i_opt];
        }
        else if (arg == "-e" && i_opt + 1 < argc)
        {
            eval_method = argv[++i_opt];
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
            colname_list_path = argv[++i_opt];
        }
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

    if (contig_list_path.empty())
    {
        std::cerr << "ERROR: contig path option is missing or failed to parse" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (eval_method.empty())
    {
        std::cerr << "ERROR: evaluate method option is missing or failed to parse" << std::endl;
        exit(EXIT_FAILURE);
    }
}

#endif //KAMRAT_EVALUATE_PARSEOPTPRINTINFO_HPP