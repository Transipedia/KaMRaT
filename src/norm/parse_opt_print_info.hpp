#ifndef KAMRAT_NORM_PARSEOPTPRINTINFO_HPP
#define KAMRAT_NORM_PARSEOPTPRINTINFO_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratNorm helper =========" << std::endl;
    std::cerr << "[Usage]      kamratQuery [-h] [-d sample_info_path] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Parameter]  -h         Print the helper" << std::endl;
    std::cerr << "             -d STRING  Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                        if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl
              << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const std::string &sample_info_path,
                         const std::string &kmer_count_path)
{
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample-info path: " << sample_info_path << std::endl;
    }
    std::cerr << "k-mer count path: " << kmer_count_path << std::endl
              << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &sample_info_path,
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
        else if (arg == "-d" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
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

    PrintRunInfo(sample_info_path, kmer_count_path);
}

#endif //KAMRAT_NORM_PARSEOPTPRINTINFO_HPP