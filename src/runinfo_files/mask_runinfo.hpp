#ifndef KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP

void MaskWelcome()
{
    std::cerr << "KaMRaT mask: k-mer sequence masking" << std::endl
              << "----------------------------------------------------------------------------------------------" << std::endl;
}

inline void PrintMaskHelper()
{
    std::cerr << "[Usage]    kamrat mask -idxdir STR -fasta STR [-options]" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -idxdir STR      Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file as the mask, mandatory" << std::endl;
    std::cerr << "            -reverse         Reverse mask, to select the k-mers in sequence fasta file [false]" << std::endl;
    std::cerr << "            -outpath STR     Path to extension results" << std::endl
              << "                                 if not provided, output to screen" << std::endl;
    std::cerr << "            -withcounts      Output sample count vectors [false]" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &idx_dir, const size_t k_len, const bool stranded,
                         const std::string &mask_file_path, const bool reverse_mask,
                         const std::string &out_path, const bool with_counts)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                  " << idx_dir << std::endl;
    std::cerr << "k-mer length:                  " << k_len << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Path to mask sequence file:    " << mask_file_path << std::endl;
    std::cerr << "Select k-mer in mask:          " << (reverse_mask ? "True" : "False") << std::endl;
    std::cerr << "Output:                        " << (out_path.empty() ? "to screen" : out_path) << ", ";
    std::cerr << (with_counts ? "with" : "without") << " count vectors" << std::endl
              << std::endl;
}

inline void ParseOptions(int argc, char *argv[],
                         std::string &idx_dir,
                         std::string &mask_file_path, bool &reverse_mask,
                         std::string &out_path, bool &with_counts)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintMaskHelper();
        exit(EXIT_SUCCESS);
    }
    std::string arg;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        arg = argv[i_opt];
        if (arg == "-h" || arg == "-help")
        {
            PrintMaskHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-idxdir" && i_opt + 1 < argc)
        {
            idx_dir = argv[++i_opt];
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            mask_file_path = argv[++i_opt];
        }
        else if (arg == "-reverse")
        {
            reverse_mask = true;
        }
        else if (arg == "-outpath" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else if (arg == "-withcounts")
        {
            with_counts = true;
        }
        else
        {
            PrintMaskHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintMaskHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (idx_dir.empty())
    {
        PrintMaskHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    if (mask_file_path.empty())
    {
        PrintMaskHelper();
        throw std::invalid_argument("-fasta STR is mandatory");
    }
}

#endif //KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP
