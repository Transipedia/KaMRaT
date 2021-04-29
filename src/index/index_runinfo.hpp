#ifndef KAMRAT_INDEX_INDEXRUNINFO_HPP
#define KAMRAT_INDEX_INDEXRUNINFO_HPP

void IndexWelcome()
{
    std::cerr << "KaMRaT index: index count table on disk" << std::endl
              << "---------------------------------------------------------------------------------" << std::endl;
}

void PrintIndexHelper()
{
    std::cerr << "[USAGE]    kamrat index -intab STR -outdir STR [-klen INT -unstrand]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]   -h, -help      Print the helper" << std::endl;
    std::cerr << "           -intab STR     Input table for index, mandatory" << std::endl;
    std::cerr << "           -outdir STR    Output index directory, mandatory" << std::endl;
    std::cerr << "           -klen          k-mer length, mandatory if features are k-mer" << std::endl
              << "                              if present, indexation will be switched to k-mer mode" << std::endl;
    std::cerr << "           -unstrand      unstranded mode, indexation with canonical k-mers" << std::endl
              << "                              if present, indexation will be switched to k-mer mode" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &count_tab_path, const std::string &out_dir,
                  const bool kmer_mode, const size_t k_len, const bool stranded)
{
    std::cerr << "Count table path:          " << count_tab_path << std::endl;
    std::cerr << "Output index directory:    " << out_dir << std::endl;
    if (kmer_mode)
    {
        std::cerr << "k-mer length:              " << k_len << std::endl;
        std::cerr << "Stranded k-mers:           " << (stranded ? "TRUE" : "FALSE") << std::endl;
    }
    else
    {
        std::cerr << "Indexation with general feature" << std::endl;
    }
    std::cerr << std::endl;
}

void ParseOptions(int argc, char *argv[], std::string &count_tab_path, std::string &out_dir,
                  bool &kmer_mode, size_t &k_len, bool &stranded)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintIndexHelper();
        exit(EXIT_SUCCESS);
    }
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
    if (kmer_mode && 0 == k_len)
    {
        PrintIndexHelper();
        throw std::invalid_argument("k-mer length in mandatory if indexing in k-mer mode");
    }
}

#endif //KAMRAT_INDEX_INDEXRUNINFO_HPP