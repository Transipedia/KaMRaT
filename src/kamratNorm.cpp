#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "run_info_parser/norm.hpp"

class SampleSum
{
public:
    SampleSum(const std::string &sample_name);
    void AddCount(double count);
    double GetCount() const;
    void Print(std::ostream &out_s) const;

private:
    std::string sample_name_;
    double sum_count_;
};

SampleSum::SampleSum(const std::string &sample_name)
    : sample_name_(sample_name), sum_count_(0)
{
}

void SampleSum::AddCount(const double count)
{
    sum_count_ += count;
}

void SampleSum::Print(std::ostream &out_s) const
{
    out_s << sample_name_ << "\t" << sum_count_ << std::endl;
}

double SampleSum::GetCount() const
{
    return sum_count_;
}

void CalcSum(std::vector<SampleSum> &sum_counts, const std::string &raw_counts_path)
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
    {
        std::istringstream conv(line);
        conv >> str_x; // feature column
        while (conv >> str_x)
        {
            sum_counts.push_back(SampleSum(str_x));
        }
    }
    std::cerr << "\tParsed sample name: " << sum_counts.size() << std::endl;

    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        size_t ns(0);
        conv >> str_x; // feature column
        while (conv >> str_x)
        {
            sum_counts[ns].AddCount(std::stod(str_x));
            ns++;
        }
    }
    raw_counts_file.close();
}

void PrintNorm(const std::vector<SampleSum> &sum_counts, const std::string &raw_counts_path, const size_t baseN)
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
    std::cout << line << std::endl;
    size_t n_line(0);
    while (std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        size_t ns(0);
        conv >> str_x;
        std::cout << str_x;
        while (conv >> str_x)
        {
            std::cout << "\t" << baseN / sum_counts.at(ns).GetCount() * std::stod(str_x);
            ns++;
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

    std::vector<SampleSum> sum_counts;
    CalcSum(sum_counts, raw_counts_path);

    PrintNorm(sum_counts, raw_counts_path, baseN);

    return EXIT_SUCCESS;
}
