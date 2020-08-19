#include <sstream>
#include <iostream>

#include "count_tab.hpp"

inline size_t WriteCountToIndex(std::ofstream &out, std::vector<float> &counts)
{
    size_t pos(out.tellp());
    out.write(reinterpret_cast<char *>(&counts[0]), counts.size() * sizeof(counts[0]));
    return pos;
}

inline void LoadCountFromIndex(std::vector<float> &counts, std::ifstream &index_file, const size_t position_on_disk, const size_t nb_samples)
{
    counts.resize(nb_samples);
    index_file.seekg(position_on_disk);
    index_file.read(reinterpret_cast<char *>(&counts[0]), nb_samples * sizeof(counts[0]));
}

inline void AddCountFromIndex(std::vector<float> &sum, std::ifstream &index_file, const size_t position_on_disk, const size_t nb_samples)
{
    std::vector<float> counts;
    LoadCountFromIndex(counts, index_file, position_on_disk, nb_samples);
    for (size_t i(0); i < nb_samples; ++i)
    {
        sum.at(i) += counts.at(i);
    }
}

CountTab::CountTab(const std::string &mode)
    : mode_(mode)
{
}

const std::string CountTab::GetMode() const
{
    return mode_;
}

const float CountTab::AddCountInMem(const std::string &line_str)
/* ------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
    Value:  inserted k-mer score (a certain column's value or -1)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. return k-mer's corresponded score value
\* ------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, line_str);
    if (nb_value_ > 0)
    {
        value_tab_.emplace_back(value_vect);
    }
    count_tab_.emplace_back(count_vect);
    return (score_pos == 0 ? -1 : value_vect[colserial_vect_[score_pos]]);
}

const float CountTab::AddIndexOnDsk(const std::string &line_str, std::ofstream &index_file)
/* -------------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
            3. indexing file stream reference
    Value:  inserted k-mer score (a certain column's value or -1)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. index the k-mer count vector in the index file and insert value vector to value table
            4. return k-mer's corresponded score value
\* -------------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, line_str);
    if (nb_value_ > 0)
    {
        value_tab_.emplace_back(value_vect);
    }
    index_pos_.emplace_back(WriteCountToIndex(index_file, count_vect));
    return (score_pos == 0 ? -1 : value_vect[colserial_vect_[score_pos]]);
}

const bool CountTab::IndexWithString(float &score_value, std::vector<float> &count_vect, const std::string &line_str, const size_t idx_pos)
/* ------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
    Value:  inserted k-mer score (a certain column's value or -1)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. return k-mer's corresponded score value
\* ------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, line_str);
    index_pos_.emplace_back(idx_pos);
    if (score_pos == 0)
    {
        score_value = 0;
        return false;
    }
    else
    {
        score_value = value_vect[colserial_vect_[score_pos]];
        return true;
    }
}

const float CountTab::GetValue(const size_t kmer_serial, const size_t valcol_serial) const
/* ------------------------------------------------------------------------------ *\
    Arg:    1. k-mer serial
            2. value column serial number
    Value:  float for the value in query
    Func:   return the value in query according to k-mer and value column serial
\* ------------------------------------------------------------------------------ */
{
    return value_tab_.at(kmer_serial).at(valcol_serial);
}

const void CountTab::GetCountInMem(std::vector<float> &count_vect, const size_t kmer_serial) const
/* ------------------------------------------------------------------ *\
    Arg:    k-mer serial
    Value:  count vector for k-mer count in query (as refArg)
    Func:   return the k-mer counts stored in memory by k-mer serial
\* ------------------------------------------------------------------ */
{
    count_vect = count_tab_.at(kmer_serial);
}

const void CountTab::GetCountOnDsk(std::vector<float> &count_vect, const size_t kmer_serial, std::ifstream &index_file) const
/* --------------------------------------------------------------------- *\
    Arg:    1. k-mer serial
            2. indexing file stream reference
    Value:  count vector for k-mer count in query (as refArg)
    Func:   return the k-mer counts stored on disk by k-mer unique code
\* --------------------------------------------------------------------- */
{
    LoadCountFromIndex(count_vect, index_file, index_pos_.at(kmer_serial), nb_count_);
}

const size_t CountTab::GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set) const
/* ------------------------------------------------------------------------------------------------------ *\
    Arg:    k-mer serial set
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. size_t for number of k-mer found in the table
    Func:   return the number of found k-mers and average counts of given k-mer list by searching memory
\* ------------------------------------------------------------------------------------------------------ */
{
    size_t nb_kmer_found(0);
    count_avg_vect.assign(nb_count_, 0);
    for (auto kmer_serial : kmer_serial_set)
    {
        ++nb_kmer_found;
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_avg_vect.at(i) += count_tab_.at(kmer_serial).at(i);
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_avg_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

const size_t CountTab::GetAvgCountOnDsk(std::vector<float> &count_avg_vect,
                                        const std::set<size_t> kmer_serial_set,
                                        std::ifstream &index_file) const
/* ------------------------------------------------------------------------------------------------------ *\
    Arg:    1. k-mer serial set
            2. indexing file stream reference
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. size_t for number of k-mer found in the table
    Func:   return the number of found k-mers and average counts of given k-mer list by searching memory
\* ------------------------------------------------------------------------------------------------------ */
{
    size_t nb_kmer_found(0);
    count_avg_vect.assign(nb_count_, 0);
    for (auto kmer_serial : kmer_serial_set)
    {
        ++nb_kmer_found;
        AddCountFromIndex(count_avg_vect, index_file, index_pos_.at(kmer_serial), nb_count_);
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_avg_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

const void CountTab::PrintFromIndex(const size_t row_serial, std::ifstream &index_file, const bool with_first_col) const
{
    size_t disk_pos = index_pos_.at(row_serial);
    std::string str_line;
    index_file.seekg(disk_pos);
    std::getline(index_file, str_line);
    if (with_first_col)
    {
        std::cout << str_line << std::endl;
    }
    else
    {
        std::cout << str_line.substr(str_line.find_first_of(" \t") + 1) << std::endl;
    }
}
