#ifndef KAMRAT_UTILS_UTILS_HPP
#define KAMRAT_UTILS_UTILS_HPP

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "sample_info.hpp"

inline void ExitIf(const bool assert_to_exit, const std::string &error_message)
{
    if (assert_to_exit)
    {
        std::cerr << error_message << std::endl;
        exit(EXIT_FAILURE);
    }
}

inline void LoadSampleInfo(SampleInfo &sample_info,
                           const std::string &sample_info_path)
{
    std::ifstream sample_condition_file(sample_info_path);
    if (!sample_condition_file.is_open())
    {
        std::cerr << "ERROR: meta data file " << sample_info_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    std::unordered_map<std::string, size_t> condition_dict;
    size_t class_num(0), label;
    while (std::getline(sample_condition_file, line))
    {
        sample_info.AddSampleInfo(line);
    }
    sample_condition_file.close();
}

inline size_t ParseHeader(std::vector<int> &header_term_nat,
                          const std::string &line,
                          const SampleInfo &sample_info,
                          const std::string &rep_value_cname)
{
    size_t nb_sample(0);
    std::istringstream conv(line);
    std::string term;
    conv >> term; // skip the first column name (sequence)
    while (conv >> term)
    {
        if (sample_info.IsEmpty() || sample_info.IsSample(term))
        {
            header_term_nat.push_back(sample_info.GetLabel(term)); // 0, 1, 2, ... for sample
            ++nb_sample;
        }
        else if (!rep_value_cname.empty() && term == rep_value_cname)
        {
            header_term_nat.push_back(-1); // -1 for rep value name
        }
        else
        {
            header_term_nat.push_back(-2); // -2 for others
        }
    }
    return nb_sample;
}

inline size_t WriteCountToIndex(std::ofstream &out, std::vector<double> &counts)
{
    size_t pos(out.tellp());
    out.write(reinterpret_cast<char *>(&counts[0]), counts.size() * sizeof(counts[0]));
    return pos;
}

inline void LoadCountFromIndex(std::ifstream &in, std::vector<double> &counts, const size_t position_on_disk, const size_t nb_samples)
{
    counts.resize(nb_samples);
    in.seekg(position_on_disk);
    in.read(reinterpret_cast<char *>(&counts[0]), nb_samples * sizeof(counts[0]));
}

inline void AddCountFromIndex(std::ifstream &in, std::vector<double> &sum, const size_t position_on_disk, const size_t nb_samples)
{
    std::vector<double> counts;
    LoadCountFromIndex(in, counts, position_on_disk, nb_samples);
    // std::cout << sum.at(0) << " += " << counts.at(0) << std::endl;
    for (size_t i(0); i < nb_samples; ++i)
    {
        sum.at(i) += counts.at(i);
    }
}

inline void GetStringLineFromDisk(std::string &str_line, std::ifstream &in_fs, const size_t disk_pos)
{
    in_fs.seekg(disk_pos);
    std::getline(in_fs, str_line);
}

#endif //KAMRAT_UTILS_UTILS_HPP