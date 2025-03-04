#ifndef KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP

const std::unordered_set<std::string> kOutFmtUniv{"tab", "bin"};
const std::unordered_set<std::string> kValueModeUniv{"int", "float"};

void MaskWelcome()
{
    std::cerr << "KaMRaT mask: mask k-mers from matrix" << std::endl
              << "----------------------------------------------------------------------------------------------" << std::endl;
}

inline void PrintMaskHelper()
{
    std::cerr << "[USAGE]    kamrat mask -idxdir STR -fasta STR [-reverse -outfmt STR -outpath STR -counts STR]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -idxdir STR      Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file as the mask, mandatory" << std::endl;
    std::cerr << "            -reverse         Reverse mask, to select the k-mers in sequence fasta file [false]" << std::endl;
    std::cerr << "            -outfmt STR      Output format, STR can be `tab` or `bin` [default `tab`]" << std::endl
              << "                                 `tab` will output the final count table, set by default" << std::endl
              << "                                 `bin` will output a binary file, to be taken by the `-with` option of other modules" << std::endl;
    std::cerr << "            -outpath STR     Path to extension results" << std::endl
              << "                                 if not provided, output to screen" << std::endl;
    std::cerr << "            -counts STR      STR can be `int` or `float` [default `int`], only works if `-outfmt tab`" << std::endl
              << "                                 `int` will round the count values to nearest integers" << std::endl
              << "                                 `float` will output the values in decimals" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &idx_dir, const size_t k_len, const bool stranded,
                         const std::string &mask_file_path, const bool reverse_mask,
                         const std::string &out_fmt, const std::string &out_path,
                         const std::string &value_mode)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                  " << idx_dir << std::endl;
    std::cerr << "k-mer length:                  " << k_len << std::endl;
    std::cerr << "Stranded mode:                 " << (stranded ? "On" : "Off") << std::endl;
    std::cerr << "Path to mask sequence file:    " << mask_file_path << std::endl;
    std::cerr << "Select k-mer in mask:          " << (reverse_mask ? "True" : "False") << std::endl;
    std::cerr << "Output:                        " << (out_path.empty() ? "to screen" : out_path) << std::endl
              << "    format:                    " + out_fmt << std::endl
              << "    value mode:                " + value_mode << std::endl
              << std::endl;
}

inline void ParseOptions(int argc, char *argv[],
                         std::string &idx_dir,
                         std::string &mask_file_path, bool &reverse_mask,
                         std::string &out_fmt, std::string &out_path,
                         std::string &value_mode)
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
        else if (arg == "-outfmt" && i_opt + 1 < argc)
        {
            out_fmt = argv[++i_opt];
        }
        else if (arg == "-outpath" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else if (arg == "-counts" && i_opt + 1 < argc)
        {
            value_mode = argv[++i_opt];
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
    if (kOutFmtUniv.find(out_fmt) == kOutFmtUniv.cend())
    {
        PrintMaskHelper();
        throw std::invalid_argument("unknown output mode: " + out_fmt);
    }
    if (kValueModeUniv.find(value_mode) == kValueModeUniv.cend())
    {
        PrintMaskHelper();
        throw std::invalid_argument("unknown value mode: " + value_mode);
    }
}

#endif // KAMRAT_RUNINFOFILES_MASKRUNINFO_HPP
