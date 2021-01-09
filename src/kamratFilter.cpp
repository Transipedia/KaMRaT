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

#include "common/tab_header.hpp"
#include "filter/filter_runinfo.hpp"

void ScanCountTab(TabHeader &tab_header,
                  const std::string &count_tab_path,
                  const size_t up_min_rec, const size_t up_min_abd,
                  const size_t down_min_rec, const size_t down_max_abd,
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

    std::string line, term;
    std::getline(kmer_count_instream, line);
    std::istringstream conv(line);
    tab_header.MakeColumns(conv, "");
    std::cout << line << std::endl;
    conv.clear();
    while (std::getline(kmer_count_instream, line))
    {
        conv.str(line);
        size_t up_rec(0), down_rec(0);
        for (size_t i(0); conv >> term && i < tab_header.GetNbCol(); ++i)
        {
            if (tab_header.IsColCount(i))
            {
                if (tab_header.GetColCondiAt(i) == "UP" && std::stof(term) >= up_min_abd)
                {
                    ++up_rec;
                }
                else if (tab_header.GetColCondiAt(i) == "DOWN" && std::stof(term) <= down_max_abd)
                {
                    ++down_rec;
                }
            }
        }
        if (up_rec >= up_min_rec && down_rec >= down_min_rec)
        {
            std::cout << line << std::endl;
        }
        conv.clear();
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
    std::string filter_info_path, count_tab_path, out_path;
    size_t up_min_rec(1), up_min_abd(1), down_min_rec(1), down_max_abd(std::numeric_limits<size_t>::max());

    ParseOptions(argc, argv, filter_info_path, up_min_rec, up_min_abd, down_min_rec, down_max_abd, out_path, count_tab_path);
    PrintRunInfo(filter_info_path, up_min_rec, up_min_abd, down_min_rec, down_max_abd, out_path, count_tab_path);

    TabHeader tab_header(filter_info_path);
    for (const auto &condi : tab_header.GetCondiNameVect())
    {
        if (condi != "UP" && condi != "DOWN")
        {
            PrintFilterHelper();
            throw std::domain_error("condition name in filter-info file could only be either UP or DOWN");
        }
    }
    ScanCountTab(tab_header, count_tab_path, up_min_rec, up_min_abd, down_min_rec, down_max_abd, out_path);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
