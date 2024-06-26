#ifndef KAMRAT_RUNINFOFILES_INDEXRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_INDEXRUNINFO_HPP

void IndexWelcome()
{
    std::cerr << "KaMRaT index: index feature count table on disk" << std::endl
              << "------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintIndexHelper()
{
    std::cerr << "[USAGE]    kamrat index -intab STR -outdir STR [-klen INT -unstrand -nfbase INT]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]   -h, -help      Print the helper" << std::endl;
    std::cerr << "           -intab STR     Input table for index, mandatory" << std::endl;
    std::cerr << "           -outdir STR    Output index directory, mandatory" << std::endl;
    std::cerr << "           -klen          k-mer length, mandatory if features are k-mer" << std::endl
              << "                              if present, indexation will be switched to k-mer mode" << std::endl
              << "                              currently, KaMRaT only supports k-mers no longer than 31nt" << std::endl
              << "                              we suggest to select k as an odd number" << std::endl;
    std::cerr << "           -unstrand      Unstranded mode, indexation with canonical k-mers" << std::endl
              << "                              if present, indexation will be switched to k-mer mode" << std::endl;
    std::cerr << "           -nfbase INT    Base for calculating normalization factor, not compatible with -nffile STR" << std::endl
              << "                              normCount_ij <- INT * rawCount_ij / sum_i{rawCount_ij}" << std::endl
              << "                              if not provided, input counts will not be normalized" << std::endl
              << "           -nffile STR    File for loading normalization factor, not compatible with -nfbase INT" << std::endl
              << "                              a tab-separated row of normalization factors, same order as table header" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &count_tab_path, const std::string &out_dir,
                  const size_t k_len, const bool stranded,
                  const size_t nf_base, const std::string &nf_file_path)
{
    std::cerr << "Count table path:             " << count_tab_path << std::endl;
    std::cerr << "Output index directory:       " << out_dir << std::endl;
    if (k_len > 0)
    {
        std::cerr << "k-mer length:                 " << k_len << std::endl;
        std::cerr << "Stranded k-mers:              " << (stranded ? "True" : "False") << std::endl;
    }
    else
    {
        std::cerr << "Indexation with general feature" << std::endl;
    }
    if (nf_base > 0)
    {
        std::cerr << "Normalization base:           " << nf_base << std::endl;
    }
    else if (!nf_file_path.empty())
    {
        std::cerr << "Normalization factor from:    " << nf_file_path << std::endl;
    }
    std::cerr << std::endl;
}

void ParseOptions(int argc, char *argv[], std::string &count_tab_path, std::string &out_dir,
                  size_t &k_len, bool &stranded, size_t &nf_base, std::string &nf_file_path)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintIndexHelper();
        exit(EXIT_SUCCESS);
    }
    bool kmer_mode(false);
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintIndexHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-intab" && i_opt + 1 < argc)
        {
            count_tab_path = argv[++i_opt];
        }
        else if (arg == "-outdir" && i_opt + 1 < argc)
        {
            out_dir = argv[++i_opt];
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_len = std::stoul(argv[++i_opt]);
            kmer_mode = true;
        }
        else if (arg == "-unstrand")
        {
            stranded = false;
            kmer_mode = true;
        }
        else if (arg == "-nfbase" && i_opt + 1 < argc)
        {
            if (!nf_file_path.empty()) {
                PrintIndexHelper();
                throw std::invalid_argument("-nfbase and -nffile cannot be given together");
            }
            nf_base = std::stoul(argv[++i_opt]);
        }
        else if (arg == "-nffile" && i_opt + 1 < argc)
        {
            if (nf_base > 0) {
                PrintIndexHelper();
                throw std::invalid_argument("-nfbase and -nffile cannot be given together");
            }
            nf_file_path = argv[++i_opt];
        }
        else
        {
            PrintIndexHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintIndexHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (count_tab_path.empty())
    {
        PrintIndexHelper();
        throw std::invalid_argument("input table is mandatory");
    }
    if (out_dir.empty())
    {
        PrintIndexHelper();
        throw std::invalid_argument("output directory for index is mandatory");
    }
    if (k_len >= 32)
    {
        PrintIndexHelper();
        throw std::invalid_argument("current KaMRaT only supports k-mer length <= 31");
    }
    if (kmer_mode && 0 == k_len)
    {
        PrintIndexHelper();
        throw std::invalid_argument("k-mer length in mandatory if indexing in k-mer mode");
    }
}

#endif //KAMRAT_RUNINFOFILES_INDEXRUNINFO_HPP