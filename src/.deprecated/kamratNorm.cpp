#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "run_info_parser/norm.hpp"
#include "data_struct/tab_header.hpp"
#include "data_struct/sample_info.hpp"

void CalcSum(TabHeader &tab_header, smpinfovect_t &smp_info_vect, const std::string &raw_counts_path)
{
    std::ifstream raw_counts_file(raw_counts_path);
    if (!raw_counts_file.is_open())
    {
        throw std::domain_error("count table " + raw_counts_path + " was not found");
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = raw_counts_path.find_last_of(".");
        if (pos != std::string::npos && raw_counts_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(raw_counts_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line, str_x;
    std::getline(kmer_count_instream, line);
    tab_header.MakeColumnInfo(line, "");
    for (size_t i(0); i < tab_header.GetNbCol(); ++i)
    {
        if (tab_header.IsSample(i))
        {
            smp_info_vect.push_back(SampleInfo(tab_header.GetColNameAt(i)));
        }
    }
    std::cerr << "\t => Number of sample parsed: " << smp_info_vect.size() << std::endl;

    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        for (size_t i(0); conv >> str_x && i < tab_header.GetNbCol(); ++i)
        {
            if (tab_header.IsSample(i))
            {
                smp_info_vect.at(tab_header.GetColSerialAt(i)).AddCount(std::stod(str_x));
            }
            if (conv.fail() || i >= tab_header.GetNbCol())
            {
                throw std::domain_error("count line string not coherent with header");
            }
        }
    }
    raw_counts_file.close();
}

void PrintSum(const std::string &smp_sum_outpath,
              const smpinfovect_t &smp_info_vect)
{
    std::ofstream sum_out(smp_sum_outpath);
    for (const auto &elem : smp_info_vect)
    {
        sum_out << elem.GetName() << "\t" << elem.GetCount() << std::endl;
    }
    sum_out.close();
}

void PrintNorm(smpinfovect_t &smp_info_vect,
               const TabHeader &tab_header,
               const std::string &raw_counts_path,
               const size_t baseN,
               const bool ln_trans,
               const std::string &out_path)
{
    std::ifstream raw_counts_file(raw_counts_path);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = raw_counts_path.find_last_of(".");
        if (pos != std::string::npos && raw_counts_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(raw_counts_file);
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

    // Parsing header line to get sample number //
    std::string line, str_x;
    std::getline(kmer_count_instream, line);
    std::cout << line << std::endl; // header line
    size_t n_line(0);
    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        for (int i(0); conv >> str_x; ++i)
        {
            if (i > 0)
            {
                std::cout << "\t";
            }
            if (tab_header.IsSample(i))
            {
                double sum_count = smp_info_vect.at(tab_header.GetColSerialAt(i)).GetCount(),
                       out_count = (sum_count == 0) ? 0 : (baseN / sum_count * std::stod(str_x));
                if (ln_trans)
                {
                    out_count = log(out_count + 1);
                }
                std::cout << out_count;
            }
            else
            {
                std::cout << str_x;
            }
        }
        std::cout << std::endl;
        n_line++;
    }

    raw_counts_file.close();
    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
}

int NormMain(int argc, char **argv)
{
    size_t baseN(0);
    std::string smp_info_path, raw_counts_path, smp_sum_path("./sample_sum.tsv"), out_path;
    bool ln_trans(false);

    ParseOptions(argc, argv, baseN, smp_info_path, ln_trans, smp_sum_path, out_path, raw_counts_path);
    PrintRunInfo(baseN, smp_info_path, ln_trans, smp_sum_path, out_path, raw_counts_path);

    smpinfovect_t smp_info_vect;
    TabHeader tab_header(smp_info_path);
    CalcSum(tab_header, smp_info_vect, raw_counts_path);
    PrintSum(smp_sum_path, smp_info_vect);

    PrintNorm(smp_info_vect, tab_header, raw_counts_path, baseN, ln_trans, out_path);

    return EXIT_SUCCESS;
}
