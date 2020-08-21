#ifndef KAMRAT_RUNINFOPARSER_MASK_HPP
#define KAMRAT_RUNINFOPARSER_MASK_HPP

#include <iostream>
#include <string>

inline void PrintHelper()
{
    std::cerr << "========= kamratMask helper =========" << std::endl;
    std::cerr << "[Usage]    kamratMask -klen INT -fasta STR [-unstranded] [-reverse-mask] KMER_COUNT_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -klen INT        Length of k-mers [31]" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file as the mask" << std::endl;
    std::cerr << "            -unstranded      If k-mers are generated from unstranded RNA-seq data" << std::endl;
    std::cerr << "            -reverse-mask    Reverse mask, to select the k-mers in sequence fasta file" << std::endl;
}

inline void PrintRunInfo(const size_t k_length,
                         const std::string &mask_file_path,
                         const bool stranded,
                         const bool reverse_mask)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
    std::cerr << "Mask sequence file:            " << mask_file_path << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Select k-mer in mask:          " << (reverse_mask ? "True" : "False") << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &k_length,
                         std::string &mask_file_path,
                         bool &stranded,
                         bool &reverse_mask,
                         std::string &count_tab_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h" || arg == "-help")
        {
            PrintHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            mask_file_path = argv[++i_opt];
        }
        else if (arg == "-unstranded")
        {
            stranded = false;
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-reverse-mask")
        {
            reverse_mask = true;
        }
        else
        {
            PrintHelper();
            throw std::domain_error("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintHelper();
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];
    if (k_length == 0)
    {
        PrintHelper();
        throw std::invalid_argument("k-mer length is missing");
    }
    if (mask_file_path.empty())
    {
        PrintHelper();
        throw std::domain_error("Mask sequence fasta path is mandatory");
    }
}

#endif //KAMRAT_RUNINFOPARSER_MASK_HPP