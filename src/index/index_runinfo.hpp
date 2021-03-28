#ifndef KAMRAT_INDEX_INDEXRUNINFO_HPP
#define KAMRAT_INDEX_INDEXRUNINFO_HPP

const void IndexWelcome()
{
    std::cerr << "KaMRaT index: index count table on disk" << std::endl
              << "----------------------------------------------------------" << std::endl;
}

const void PrintIndexHelper()
{
    std::cerr << "[USAGE]    kamrat index -out-dir STR INPUT_TAB_PATH" << std::endl
              << std::endl;
    std::cerr << "[OPTION]   -h, -help    Print the helper" << std::endl;
    std::cerr << "           -out-dir     Output directory for storing index" << std::endl
              << std::endl;
}

const void PrintRunInfo(const std::string &out_dir, const std::string &count_tab_path)
{
    std::cerr << "Output directory for index:    " << std::endl;
    std::cerr << "Count table path:              " << count_tab_path << std::endl;
    std::cerr << std::endl;
}

const void ParseOptions(int argc, char *argv[],
                        std::string &out_dir,
                        std::string &count_tab_path)
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
        else if (arg == "-out-dir" && i_opt + 1 < argc)
        {
            out_dir = argv[++i_opt];
        }
        else
        {
            PrintIndexHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintIndexHelper();
        throw std::domain_error("k-mer count table path is mandatory");
    }
    count_tab_path = argv[i_opt++];
    if (out_dir.empty())
    {
        PrintIndexHelper();
        throw std::domain_error("output directory for index is mandatory");
    }
}

#endif //KAMRAT_INDEX_INDEXRUNINFO_HPP