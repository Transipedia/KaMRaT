#ifndef KAMRAT_QUERY_PARSEOPTPRINTINFO_HPP
#define KAMRAT_QUERY_PARSEOPTPRINTINFO_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratQuery helper =========" << std::endl;
    std::cerr << "[Usage]    kamratQuery -l contig_list_path [-n] [-k k_length] [-d colname_list_path] [-m query_mode] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h         Print the helper" << std::endl;
    std::cerr << "            -l STRING  Contig list file path (MANDATORY)" << std::endl
              << "                       could be fasta, list, or table with contig sequences as the first column" << std::endl;
    std::cerr << "            -k INT     k-mer length [31]" << std::endl;
    std::cerr << "            -n         If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -d STRING  Colname list path indicating columns to print, either list or table with sample names as the first column" << std::endl
              << "                       if absent, all columns will be printed" << std::endl;
    std::cerr << "            -m STRING  Query mode (extract, eliminate, or estimate), [extract]" << std::endl
              << "                       extract for extracting from the k-mer count table thoses k-mers that are present in contigs" << std::endl
              << "                       eliminate for eliminating from the k-mer count table thoses k-mers that are present in contigs" << std::endl
              << "                       estimate for estimating statistics (mean, sd, etc.) from the k-mer count table contig by contig" << std::endl
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

inline void PrintRunInfo(const std::string &contig_list_path,
                         const int k_length,
                         const bool stranded,
                         const std::string &colname_list_path,
                         const std::string &query_mode,
                         const std::string &estimate_method,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "Contig list path:       " << contig_list_path << std::endl;
    std::cerr << "Stranded mode:          " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:           " << k_length << std::endl;
    if (!colname_list_path.empty())
    {
        std::cerr << "Colname list path:      " << colname_list_path << std::endl;
    }
    std::cerr << "Query mode:             " << query_mode << std::endl;
    if (query_mode == "estimate")
    {
        std::cerr << "Estimate method:        " << estimate_method << std::endl;
    }
    std::cerr << "k-mer count path:       " << kmer_count_path << std::endl;
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &contig_list_path,
                         int &k_length,
                         bool &stranded,
                         std::string &colname_list_path,
                         std::string &query_mode,
                         std::string &estimate_method,
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
        else if (arg == "-m" && i_opt + 1 < argc)
        {
            query_mode = argv[++i_opt];
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

    SubCommandParser(query_mode, estimate_method);

    if (contig_list_path.empty())
    {
        std::cerr << "ERROR: contig path option is missing or failed to parse" << std::endl;
        exit(EXIT_FAILURE);
    }
}

#endif //KAMRAT_QUERY_PARSEOPTPRINTINFO_HPP