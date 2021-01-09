#ifndef KAMRAT_MASK_MASKRUNINFO_HPP
#define KAMRAT_MASK_MASKRUNINFO_HPP

#include <iostream>
#include <string>

inline void PrintMaskHelper()
{
    std::cerr << "[Usage]    kamrat mask -klen INT -fasta STR [-options] KMER_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -klen INT        Length of k-mers, mandatory" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file as the mask, mandatory" << std::endl;
    std::cerr << "            -unstrand        If k-mers are generated from unstranded RNA-seq data" << std::endl;
    std::cerr << "            -reverse-mask    Reverse mask, to select the k-mers in sequence fasta file" << std::endl;
    std::cerr << "            -out-path        Output table path [default: output to screen]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const size_t k_length,
                         const std::string &mask_file_path,
                         const bool stranded,
                         const bool reverse_mask,
                         const std::string &out_path)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:                  " << k_length << std::endl;
    std::cerr << "Mask sequence file:            " << mask_file_path << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Select k-mer in mask:          " << (reverse_mask ? "True" : "False") << std::endl;
    if (!out_path.empty())
    {
        std::cerr << "Output path:                   " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &k_length,
                         std::string &mask_file_path,
                         bool &stranded,
                         bool &reverse_mask,
                         std::string &out_path,
                         std::string &count_tab_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h" || arg == "-help")
        {
            PrintMaskHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_length = atoi(argv[++i_opt]);
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            mask_file_path = argv[++i_opt];
        }
        else if (arg == "-unstrand")
        {
            stranded = false;
        }
        else if (arg == "-reverse-mask")
        {
            reverse_mask = true;
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else
        {
            PrintMaskHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintMaskHelper();
        throw std::invalid_argument("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];
    if (k_length == 0)
    {
        PrintMaskHelper();
        throw std::invalid_argument("k-mer length is mandatory");
    }
    if (mask_file_path.empty())
    {
        PrintMaskHelper();
        throw std::invalid_argument("Mask sequence fasta path is mandatory");
    }
}

#endif //KAMRAT_MASK_MASKRUNINFO_HPP
