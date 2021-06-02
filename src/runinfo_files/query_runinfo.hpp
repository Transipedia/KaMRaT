#ifndef KAMRAT_RUNINFOFILES_QUERY_HPP
#define KAMRAT_RUNINFOFILES_QUERY_HPP

#include <unordered_set>

const std::unordered_set<std::string> kQueryMethodUniv{"mean", "median"};

void QueryWelcome()
{
    std::cerr << "KaMRaT query: query sequences" << std::endl
              << "-------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintQueryHelper()
{
    std::cerr << "[USAGE]    kamrat query -fasta STR -toquery STR [-options]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help         Print the helper" << std::endl;
    std::cerr << "            -idxdir STR      Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -fasta STR       Sequence fasta file, mandatory" << std::endl;
    std::cerr << "            -toquery STR     Query method, mandatory, can be one of:" << std::endl
              << "                                 mean        mean count among all composite k-mers for each sample" << std::endl
              << "                                 median      median count among all composite k-mers for each sample" << std::endl;
    std::cerr << "            -withabsent      Output also absent queries (count vector all 0) [default: false]" << std::endl;
    std::cerr << "            -outpath STR     Path to extension results" << std::endl
              << "                                 if not provided, output to screen" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir, const std::string &seq_file_path,
                  const std::string &query_mtd, const bool with_absent, const std::string &out_path)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:             " << idx_dir << std::endl;
    std::cerr << "Path to sequence file:    " << seq_file_path << std::endl;
    std::cerr << "Query method:             " << query_mtd << std::endl;
    std::cout << "Output absent query:      " << (with_absent ? "True" : "False") << std::endl;
    std::cerr << "Output:                   " << (out_path.empty() ? "to screen" : out_path) << std::endl
              << std::endl;
}

void ParseOptions(int argc, char *argv[], std::string &idx_dir, std::string &seq_file_path,
                  std::string &query_mtd, bool &with_absent, std::string &out_path)
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
}

#endif //KAMRAT_RUNINFOFILES_QUERY_HPP
