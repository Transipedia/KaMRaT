#ifndef KAMRAT_RUNINFOFILES_FILTERRUNINFO_HPP
#define KAMRAT_RUNINFOFILES_FILTERRUNINFO_HPP

const std::unordered_set<std::string> kOutFmtUniv{"tab", "fa", "bin"};
const std::unordered_set<std::string> kValueModeUniv{"int", "float"};

void FilterWelcome()
{
    std::cerr << "KaMRaT filter: filter feature by expression level" << std::endl
              << "-----------------------------------------------------------------------------------------------------------" << std::endl;
}

void PrintFilterHelper()
{
    std::cerr << "[USAGE]    kamrat filter -idxdir STR -design STR [-upmin INT1:INT2 -downmax INT1:INT2 -reverse -outfmt STR -outpath STR -counts STR]" << std::endl
              << std::endl;
    std::cerr << "[OPTION]    -h,-help              Print the helper" << std::endl;
    std::cerr << "            -idxdir STR           Indexing folder by KaMRaT index, mandatory" << std::endl;
    std::cerr << "            -design STR           Path to filter design file, a table of two columns, mandatory" << std::endl
              << "                                      the first column indicate sample names" << std::endl
              << "                                      the second column should be either UP or DOWN (capital letters)" << std::endl
              << "                                          samples with UP will be considered as up-regulated samples" << std::endl
              << "                                          samples with DOWN will be considered as down-regulated samples" << std::endl
              << "                                          samples not given will be neutral (not considered for filter)" << std::endl
              << "                                          samples can also be all UP or all DOWN" << std::endl;
    std::cerr << "            -upmin INT1:INT2      Up feature lower bound, [1:1, meaning no filter]" << std::endl
              << "                                      output features counting >= INT1 in >= INT2 UP-samples" << std::endl;
    std::cerr << "            -downmax INT1:INT2    Down feature upper bound [inf:1, meaning no filter]" << std::endl
              << "                                      output features counting <= INT1 in >= INT2 DOWN-samples" << std::endl;
    std::cerr << "            -reverse              Reverse filter, to remove eligible features [false]" << std::endl;
    std::cerr << "            -outfmt STR           Output format, STR can be `tab`, `fa`, or `bin` [default `tab`]" << std::endl
              << "                                      `tab` will output the final count table, set by default" << std::endl
              << "                                      `fa` will output a fasta file containing sequences without counts" << std::endl
              << "                                      `bin` will output a binary file, to be taken by the `-with` option of other modules" << std::endl;
    std::cerr << "            -outpath STR          Path to results after filter" << std::endl
              << "                                      if not provided, output to screen" << std::endl;
    std::cerr << "            -counts STR           STR can be `int` or `float` [default `int`], only works if `-outfmt tab`" << std::endl
              << "                                      `int` will round the count values to nearest integers" << std::endl
              << "                                      `float` will output the values in decimals" << std::endl
              << std::endl;
}

void PrintRunInfo(const std::string &idx_dir,
                  const std::string &dsgn_path,
                  const size_t up_min_abd, const size_t up_min_rec,
                  const size_t down_max_abd, const size_t down_min_rec,
                  const bool reverse_filter,
                  const std::string &out_fmt, const std::string &out_path,
                  const std::string &value_mode)
{
    std::cerr << std::endl;
    std::cerr << "KaMRaT index:                  " << idx_dir << std::endl;
    std::cerr << "Path to filter design file:    " << dsgn_path << std::endl;
    std::cerr << "Up-regulated lower bound:      " << std::endl
              << "\tfeatures counting >= " << up_min_abd << " in >= " << up_min_rec << " up-regulated samples" << std::endl;
    std::cerr << "Down-regulated upper bound:    " << std::endl
              << "\tfeatures counting <= " << (down_max_abd == std::numeric_limits<size_t>::max() ? "inf" : std::to_string(down_max_abd))
              << " in >= " << down_min_rec << " down-regulated samples" << std::endl;
    std::cerr << "Remove eligible features:      " << (reverse_filter ? "True" : "False") << std::endl;
    std::cerr << "Output:                        " << (out_path.empty() ? "to screen" : out_path) << std::endl
              << "    format:                    " + out_fmt << std::endl
              << "    value mode:                " + value_mode << std::endl
              << std::endl;
}

void ParseOptions(int argc, char *argv[],
                  std::string &idx_dir,
                  std::string &dsgn_path,
                  size_t &up_min_abd, size_t &up_min_rec,
                  size_t &down_max_abd, size_t &down_min_rec,
                  bool &reverse_filter,
                  std::string &out_fmt, std::string &out_path,
                  std::string &value_mode)
{
    int i_opt(1);
    if (argc == 1)
    {
        PrintFilterHelper();
        exit(EXIT_SUCCESS);
    }
    size_t split_pos;
    std::string arg;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        arg = argv[i_opt];
        if (arg == "-help" || arg == "-h")
        {
            PrintFilterHelper();
            exit(EXIT_SUCCESS);
        }
        if (arg == "-idxdir" && i_opt + 1 < argc)
        {
            idx_dir = argv[++i_opt];
        }
        else if (arg == "-design" && i_opt + 1 < argc)
        {
            dsgn_path = argv[++i_opt];
        }
        else if (arg == "-upmin" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                up_min_abd = std::stoul(arg.substr(0, split_pos));
                up_min_rec = std::stoul(arg.substr(split_pos + 1));
            }
            else
            {
                throw std::invalid_argument("unable to parse -upmin argument: " + arg);
            }
        }
        else if (arg == "-downmax" && i_opt + 1 < argc)
        {
            arg = argv[++i_opt];
            split_pos = arg.find(":");
            if (split_pos != std::string::npos)
            {
                down_max_abd = std::stoul(arg.substr(0, split_pos));
                down_min_rec = std::stoul(arg.substr(split_pos + 1));
            }
            else
            {
                throw std::invalid_argument("unable to parse -downmax argument: " + arg);
            }
        }
        else if (arg == "-reverse")
        {
            reverse_filter = true;
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
            PrintFilterHelper();
            throw std::invalid_argument("unknown option: " + arg);
        }
        ++i_opt;
    }
    if (i_opt < argc)
    {
        PrintFilterHelper();
        throw std::invalid_argument("cannot parse arguments after " + std::string(argv[i_opt]));
    }
    if (idx_dir.empty())
    {
        PrintFilterHelper();
        throw std::invalid_argument("-idxdir STR is mandatory");
    }
    if (dsgn_path.empty())
    {
        PrintFilterHelper();
        throw std::domain_error("-design STR is mandatory");
    }
    if (kOutFmtUniv.find(out_fmt) == kOutFmtUniv.cend())
    {
        PrintFilterHelper();
        throw std::invalid_argument("unknown output mode: " + out_fmt);
    }
    if (kValueModeUniv.find(value_mode) == kValueModeUniv.cend())
    {
        PrintFilterHelper();
        throw std::invalid_argument("unknown value mode: " + value_mode);
    }
}

#endif //KAMRAT_RUNINFOFILES_FILTERRUNINFO_HPP
