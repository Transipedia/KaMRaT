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
#include "data_struct/count_tab_header.hpp"
#include "data_struct/sample_info.hpp"

void CalcSum(sampleInfoVect_t &sample_info_vect, const std::string &sample_info_path, const std::string &raw_counts_path)
{
    std::ifstream raw_counts_file(raw_counts_path);
    if (!raw_counts_file.is_open())
    {
        throw std::domain_error("count table " + raw_counts_path + " was not found");
    }

    size_t pos = raw_counts_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && raw_counts_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(raw_counts_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line, str_x;
    std::getline(kmer_count_instream, line);
    CountTabHeader count_tab_header;
    count_tab_header.MakeColumnInfo(line, sample_info_path, "NULLFORNOSCORE");
    for (int i(0); i < count_tab_header.GetNbColumn(); ++i)
    {
        if (count_tab_header.GetColNature(i) == 's')
        {
            sample_info_vect.push_back(SampleInfo(count_tab_header.GetColName(i)));
        }
    }
    std::cerr << "\t => Number of sample parsed: " << sample_info_vect.size() << std::endl;

    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        for (int i(0); conv >> str_x; ++i)
        {
            if (count_tab_header.GetColNature(i) == 's')
            {
                sample_info_vect.at(count_tab_header.GetColSerial(i)).AddCount(std::stod(str_x));
            }
        }
    }
    raw_counts_file.close();
}

void PrintNorm(sampleInfoVect_t &sample_info_vect,
               const std::string &sample_info_path,
               const std::string &raw_counts_path,
               const size_t baseN,
               const std::string &trans_mode)
{
    std::ifstream raw_counts_file(raw_counts_path);

    size_t pos = raw_counts_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && raw_counts_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(raw_counts_file);
    std::istream kmer_count_instream(&inbuf);

    // Parsing header line to get sample number //
    std::string line, str_x;
    std::getline(kmer_count_instream, line);
    std::cout << line << std::endl; // header line
    CountTabHeader count_tab_header;
    count_tab_header.MakeColumnInfo(line, sample_info_path, "NULLFORNOSCORE");

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
            if (count_tab_header.GetColNature(i) == 's')
            {
                double sum_count = sample_info_vect.at(count_tab_header.GetColSerial(i)).GetCount(),
                       out_count = (sum_count == 0) ? 0 : (baseN / sum_count * std::stod(str_x));
                if (trans_mode == "log")
                {
                    out_count = (out_count == 0) ? 0 : log(out_count);
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

    if (raw_counts_file.is_open())
    {
        raw_counts_file.close();
    }
}

int main(int argc, char **argv)
{
    size_t baseN(0);
    std::string sample_info_path, trans_mode, raw_counts_path;

    ParseOptions(argc, argv, baseN, sample_info_path, trans_mode, raw_counts_path);
    PrintRunInfo(baseN, sample_info_path, trans_mode, raw_counts_path);

    sampleInfoVect_t sample_info_vect;
    CalcSum(sample_info_vect, sample_info_path, raw_counts_path);

    PrintNorm(sample_info_vect, sample_info_path, raw_counts_path, baseN, trans_mode);

    return EXIT_SUCCESS;
}
