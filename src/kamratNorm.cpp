#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "norm/parse_opt_print_info.hpp"
#include "norm/sample_sum.hpp"

#define NORM_BASE 1E9

void CalcSum(std::vector<SampleSum> &sum_counts, const std::string &raw_counts_path)
{
    std::ifstream raw_counts_file(raw_counts_path);
    if (!raw_counts_file.is_open())
    {
        std::cerr << "ERROR: k-mer count file " << raw_counts_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line, str_x;
    std::getline(raw_counts_file, line);
    {
        std::istringstream conv(line);
        conv >> str_x; // feature column
        while (conv >> str_x)
        {
            sum_counts.push_back(SampleSum(str_x));
        }
    }
    std::cerr << "Parsed sample name: " << sum_counts.size() << std::endl;

    while (std::getline(raw_counts_file, line))
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

    if (raw_counts_file.is_open())
    {
        raw_counts_file.close();
    }
}

void PrintNorm(const std::vector<SampleSum> &sum_counts, const std::string &raw_counts_path)
{
    std::ifstream raw_counts_file(raw_counts_path);
    size_t pos = raw_counts_path.find_last_of(".");

    // Parsing header line to get sample number //
    std::string line, str_x;
    std::getline(raw_counts_file, line);
    std::cout << line << std::endl;
    size_t n_line(0);
    while (std::getline(raw_counts_file, line))
    {
        std::istringstream conv(line);
        size_t ns(0);
        conv >> str_x;
        std::cout << str_x;
        while (conv >> str_x)
        {
            std::cout << "\t" << NORM_BASE / sum_counts[ns].GetCount() * std::stod(str_x);
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

int main(int argc, char *argv[])
{
    std::string sample_info_path, kmer_count_path;
    ParseOptions(argc, argv, sample_info_path, kmer_count_path);

    if (!sample_info_path.empty())
    {
    }

    std::ofstream fout("sum_counts.tsv");
    std::vector<SampleSum> sum_counts;
    CalcSum(sum_counts, raw_counts_path);
    for (const auto &x : sum_counts)
    {
        x.Print(fout);
    }
    fout.close();
    PrintNorm(sum_counts, raw_counts_path);

    return EXIT_SUCCESS;
}