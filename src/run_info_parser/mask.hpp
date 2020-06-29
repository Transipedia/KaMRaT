#ifndef KAMRAT_RUNINFOPARSER_MASK_HPP
#define KAMRAT_RUNINFOPARSER_MASK_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratMask helper =========" << std::endl;
    std::cerr << "[Usage]    kamratMask -f mask.seq.fasta [-n] [-k k_length] [-r] kmer.count.matrix" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h           Print the helper" << std::endl;
    std::cerr << "            -f STRING    Mask sequence fasta file" << std::endl;
    std::cerr << "            -n           If k-mers are generated from unstranded RNA-seq data" << std::endl;
    std::cerr << "            -k INT       Length of k-mers [31]" << std::endl;
    std::cerr << "            -r           Reverse mask, to select the k-mers exist in mask sequence fasta file" << std::endl;
    exit(EXIT_SUCCESS);
}

inline void PrintRunInfo(const std::string &mask_file_path,
                         const bool stranded,
                         const size_t k_length,
                         const bool reverse_mask)
{
    std::cerr << std::endl;
    std::cerr << "Mask sequence file:            " << mask_file_path << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
    std::cerr << "Select k-mer in mask:          " << (reverse_mask ? "True" : "Off") << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &mask_file_path,
                         bool &stranded,
                         size_t &k_length,
                         bool &reverse_mask,
                         std::string &count_tab_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h")
        {
            PrintHelper();
        }
        else if (arg == "-f" && i_opt + 1 < argc)
        {
            mask_file_path = argv[++i_opt];
        }
        else if (arg == "-n")
        {
            stranded = false;
        }
        else if (arg == "-k" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-r")
        {
            reverse_mask = true;
        }
        else
        {
            throw std::domain_error("unknown option " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];
    if (mask_file_path.empty())
    {
        throw std::domain_error("Mask sequence fasta path is mandatory");
    }
}

#endif //KAMRAT_RUNINFOPARSER_MASK_HPP