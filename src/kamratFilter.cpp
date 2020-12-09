#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <ctime>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "data_struct/tab_header.hpp"
#include "run_info_parser/filter.hpp"

const long GetSampleSerial(const std::string &sample_name, const TabHeader &tab_header)
{
    for (size_t i(0); i < tab_header.GetNbCol(); ++i)
    {
        if (tab_header.GetColNameAt(i) == sample_name)
        {
            return tab_header.GetColSerialAt(i);
        }
    }
    throw std::domain_error("unknown sample names for expressed sample filter: " + sample_name);
}

const bool IsSampleToCount(const std::string &this_level,
                           const std::string &this_name,
                           const char this_label,
                           const std::string &other_level,
                           const std::string &other_name,
                           const char other_label,
                           const TabHeader &tab_header,
                           const size_t i_col)
{
    if (this_name == "all") // all samples
    {
        return true;
    }
    else if (this_name == "rest" && other_level == "smp" && tab_header.GetColNameAt(i_col) != other_name) // samples other than silent sample
    {
        return true;
    }
    else if (this_name == "rest" && other_level == "cond" && tab_header.GetColNatureAt(i_col) != other_label) // simples other than silent condition
    {
        return true;
    }
    else if (this_level == "smp" && tab_header.GetColNameAt(i_col) == this_name) // expressed sample
    {
        return true;
    }
    else if (this_level == "cond" && tab_header.GetColNatureAt(i_col) == this_label) // samples in expressed condition
    {
        return true;
    }
    return false;
}

void ScanCountTab(TabHeader &tab_header,
                  const std::string &count_tab_path,
                  const std::string &express_level,
                  const std::string &express_name,
                  const size_t express_min_abd,
                  const size_t express_min_rec,
                  const std::string &silent_level,
                  const std::string &silent_name,
                  const size_t silent_max_abd,
                  const size_t silent_min_rec,
                  const std::string &out_path)
{
    std::ifstream count_tab_file(count_tab_path);
    if (!count_tab_file.is_open())
    {
        throw std::domain_error("count table file " + count_tab_path + " is not found");
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = count_tab_path.find_last_of(".");
        if (pos != std::string::npos && count_tab_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(count_tab_file);
    std::istream kmer_count_instream(&inbuf);

    std::ofstream out_file;
    if (!out_path.empty())
    {
        out_file.open(out_path);
        if (!out_file.is_open())
        {
            throw std::domain_error("cannot open file: " + out_path);
        }
    }
    auto backup_buf = std::cout.rdbuf();
    if (!out_path.empty()) // output to file if a path is given, to screen if not
    {
        std::cout.rdbuf(out_file.rdbuf());
    }

    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::set<std::string> preserved_condition_names{"all", "rest"};
    std::string line;
    char express_label = (express_level == "cond" ? tab_header.GetConditionLabel(express_name) : '\0'),
         silent_label = (silent_level == "cond" ? tab_header.GetConditionLabel(silent_name) : '\0');
    std::getline(kmer_count_instream, line);
    tab_header.MakeColumnInfo(line, "");
    std::cout << line << std::endl;
    //----- Dealing with Following k-mer Count Lines -----//
    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        std::string term;
        size_t express_rec(0), silent_rec(0);
        for (size_t i(0); conv >> term && i < tab_header.GetNbCol(); ++i)
        {
            if (!tab_header.IsSample(i))
            {
                continue;
            }
            if (IsSampleToCount(express_level, express_name, express_label, silent_level, silent_name, silent_label, tab_header, i) &&
                std::stof(term) >= express_min_abd)
            {
                ++express_rec;
            }
            if (IsSampleToCount(silent_level, silent_name, silent_label, express_level, express_name, express_label, tab_header, i) &&
                std::stof(term) <= silent_max_abd)
            {
                ++silent_rec;
            }
            if (conv.fail() || i == tab_header.GetNbCol())
            {
                throw std::domain_error("line parsing not coherent with header");
            }
        }
        // std::cerr << express_rec << "\t" << silent_rec << std::endl;
        if (express_rec >= express_min_rec && silent_rec >= silent_min_rec)
        {
            std::cout << line << std::endl;
        }
    }
    count_tab_file.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
}

int FilterMain(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string count_tab_path, sample_info_path, express_level, express_name, silent_level, silent_name, out_path;
    int express_min_abd(-1), express_min_rec(-1), silent_max_abd(-1), silent_min_rec(-1);

    ParseOptions(argc, argv, count_tab_path,
                 express_level, express_name, express_min_rec, express_min_abd,
                 silent_level, silent_name, silent_min_rec, silent_max_abd,
                 sample_info_path, out_path);
    PrintRunInfo(count_tab_path, express_level, express_name, express_min_rec, express_min_abd,
                 silent_level, silent_name, silent_min_rec, silent_max_abd, sample_info_path, out_path);

    TabHeader tab_header(sample_info_path, {"all", "rest"});
    ScanCountTab(tab_header, count_tab_path, express_level, express_name, express_min_abd, express_min_rec,
                 silent_level, silent_name, silent_max_abd, silent_min_rec, out_path);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
