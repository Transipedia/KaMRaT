#ifndef KAMRAT_RUNINFOFILES_QUERY_HPP
#define KAMRAT_RUNINFOFILES_QUERY_HPP

#include <unordered_set>

const std::unordered_set<std::string> kQueryMethodUniv{"mean", "median"};
const std::unordered_set<std::string> kValueModeUniv{"int", "float"};

void QueryWelcome()
{
    std::cerr << "KaMRaT query: query sequences counts" << std::endl
              << "-------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintQueryHelper()
{
    std::cerr << "[USAGE]    kamrat query -idxdir STR -fasta STR -toquery STR [-withabsent -outpath STR -counts]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -idxdir STR      Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file, mandatory" << std::endl;
    std::cerr << "            -toquery STR     Query method, mandatory, can be one of:" << std::endl
              << "                                 mean        mean count among all composite k-mers for each sample" << std::endl
              << "                                 median      median count among all composite k-mers for each sample" << std::endl;
    std::cerr << "            -withabsent      Output also absent queries (count vector all 0) [default: false]" << std::endl;
    std::cerr << "            -outpath STR     Path to extension results" << std::endl
              << "                                 if not provided, output to screen" << std::endl;
    std::cerr << "            -counts STR      STR can be `int` or `float` [default `int`], only works if `-outfmt tab`" << std::endl
              << "                                 `int` will round the count values to nearest integers" << std::endl
              << "                                 `float` will output the values in decimals" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir, const std::string &seq_file_path,
                  const std::string &query_mtd, const bool with_absent, const std::string &out_path,
                  const std::string &value_mode)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:             " << idx_dir << std::endl;
    std::cerr << "Path to sequence file:    " << seq_file_path << std::endl;
    std::cerr << "Query method:             " << query_mtd << std::endl;
    std::cout << "Output absent query:      " << (with_absent ? "True" : "False") << std::endl;
    std::cerr << "Output:                   " << (out_path.empty() ? "to screen" : out_path) << std::endl
              << "    value mode:           " + value_mode << std::endl
              << std::endl;
}

void ParseOptions(int argc, char *argv[], std::string &idx_dir, std::string &seq_file_path,
                  std::string &query_mtd, bool &with_absent, std::string &out_path,
                  std::string &value_mode)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintQueryHelper();
        exit(EXIT_SUCCESS);
    }
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintQueryHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-idxdir" && i_opt + 1 < argc)
        {
            idx_dir = argv[++i_opt];
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            seq_file_path = argv[++i_opt];
        }
        else if (arg == "-toquery" && i_opt + 1 < argc)
        {
            query_mtd = argv[++i_opt];
        }
        else if (arg == "-withabsent")
        {
            with_absent = true;
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
            PrintQueryHelper();
            throw std::invalid_argument("unknown option " + std::string(argv[i_opt]));
        }
        ++i_opt;
    }
    if (idx_dir.empty())
    {
        PrintQueryHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    if (seq_file_path.empty())
    {
        PrintQueryHelper();
        throw std::invalid_argument("-fasta STR is mandatory");
    }
    if (query_mtd.empty())
    {
        PrintQueryHelper();
        throw std::invalid_argument("-toquery STR is mandatory");
    }
    else if (kQueryMethodUniv.find(query_mtd) == kQueryMethodUniv.cend())
    {
        PrintQueryHelper();
        throw std::invalid_argument("unknown query method: " + query_mtd);
    }
    if (kValueModeUniv.find(value_mode) == kValueModeUniv.cend())
    {
        PrintQueryHelper();
        throw std::invalid_argument("unknown value mode: " + value_mode);
    }
}

#endif // KAMRAT_RUNINFOFILES_QUERY_HPP
