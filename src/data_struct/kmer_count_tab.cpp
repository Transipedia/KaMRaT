#include <sstream>

#include "kmer_count_tab.hpp"

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

inline void GetStringLineFromDisk(std::string &str_line, std::ifstream &index_file, const size_t disk_pos)
{
    index_file.seekg(disk_pos);
    std::getline(index_file, str_line);
}

inline void StrLine2ValueCountVects(std::vector<float> &value_vect,
                                    std::vector<float> &count_vect,
                                    const std::string &line_str,
                                    const ColumnInfo &col_info)
{
    std::istringstream conv(line_str);
    std::string term;
    for (unsigned int i(0); i < col_info.GetNbColumn() && conv >> term; ++i)
    {
        if (col_info.GetColumnNature(i) == 'v')
        {
            value_vect.push_back(std::stof(term));
        }
        else if (col_info.GetColumnNature(i) == 's')
        {
            count_vect.push_back(std::stof(term));
        }
    }
}

KMerCountTab::KMerCountTab(const std::string &mode)
    : mode_(mode), nb_value_(0), nb_count_(0)
{
}

const std::string KMerCountTab::GetMode() const
{
    return mode_;
}

const float KMerCountTab::AddKMerCountInMem(const uint64_t kmer_code,
                                            const std::string &line_str,
                                            const ColumnInfo &column_info,
                                            const std::string &score_colname)
/* ----------------------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. string of input k-mer table line
            3. ColumnInfo object
            4. column name for score
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value & count vector to insert is coherent with the tables
            2. parse the line string
            3. check whether the k-mer code has been alreaded inserted (avoiding duplicate input)
            4. insert the k-mer value & count vector to k-mer value & count table
            5. associate the k-mer unique code with its inserted row number in the value table
            6. return k-mer's corresponded score value
\* ----------------------------------------------------------------------------------------------- */
{
    std::vector<float> count_vect, value_vect;
    StrLine2ValueCountVects(value_vect, count_vect, line_str, column_info);
    if (value_tab_.empty() || count_tab_.empty())
    {
        nb_value_ = value_vect.size();
        nb_count_ = count_vect.size();
    }
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    if (!kmer_serial_.insert({kmer_code, value_tab_.size()}).second)
    {
        throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
    }
    value_tab_.emplace_back(value_vect);
    count_tab_.emplace_back(count_vect);
    return (score_colname.empty() ? value_tab_.size() : value_vect.at(column_info.GetColumnSerial(score_colname)));
}

const float KMerCountTab::AddKMerIndexOnDsk(const uint64_t kmer_code,
                                            const std::string &line_str,
                                            const ColumnInfo &column_info,
                                            const std::string &score_colname,
                                            std::ofstream &index_file)
/* -------------------------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. string of input k-mer table line
            3. ColumnInfo object
            4. column name for score
            5. indexing file stream reference
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value vector to insert is coherent with the table
            2. parse the line string
            3. check whether the k-mer code has been alreaded inserted (avoiding duplicate input)
            4. index the k-mer count vector in the index file and insert value vector to value table
            5. associate the k-mer unique code with its indexed position in the file
            6. return k-mer's corresponded score value
\* -------------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    StrLine2ValueCountVects(value_vect, count_vect, line_str, column_info);
    if (value_tab_.empty() || index_pos_.empty())
    {
        nb_value_ = value_vect.size();
        nb_count_ = count_vect.size();
    }
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    if (!kmer_serial_.insert({kmer_code, value_tab_.size()}).second)
    {
        throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
    }
    value_tab_.emplace_back(value_vect);
    index_pos_.emplace_back(WriteCountToIndex(index_file, count_vect));
    return (score_colname.empty() ? value_tab_.size() : value_vect.at(column_info.GetColumnSerial(score_colname)));
}

const bool KMerCountTab::IsKMerExist(const uint64_t kmer_code) const
/* ---------------------------------------------------------------------------- *\
    Arg:    k-mer unique code
    Value   bool indicating whether the k-mer exists in k-mer value dictionary
    Func:   check if the k-mer has been registered in the table
\* ---------------------------------------------------------------------------- */
{
    return (kmer_serial_.find(kmer_code) != kmer_serial_.cend());
}

const float KMerCountTab::GetValue(const uint64_t kmer_code, const size_t i_valcol) const
/* --------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. value name string
            3. value column serial number
    Value:  1. float for the value in query (as refArg)
            2. bool for whether the query was succeeded
    Func:   return the value in query according to k-mer unique code and value name
\* --------------------------------------------------------------------------------- */
{
    auto iter_kmer = kmer_serial_.find(kmer_code);
    if (iter_kmer == kmer_serial_.cend())
    {
        throw std::domain_error("non-registered k-mer code" + std::to_string(kmer_code));
    }
    return value_tab_.at(iter_kmer->second).at(i_valcol);
}

const bool KMerCountTab::GetCountInMem(std::vector<float> &count_vect, const uint64_t kmer_code) const
/* ----------------------------------------------------------------------- *\
    Arg:    k-mer unique code
    Value:  1. count vector for k-mer count in query (as refArg)
            2. bool for whether the query was succeeded
    Func:   return the k-mer counts stored in memory by k-mer unique code
\* ----------------------------------------------------------------------- */
{
    auto iter = kmer_serial_.find(kmer_code);
    if (iter == kmer_serial_.cend())
    {
        return false;
    }
    count_vect = count_tab_.at(iter->second);
    return true;
}

const bool KMerCountTab::GetCountOnDsk(std::vector<float> &count_vect, const uint64_t kmer_code, std::ifstream &index_file, const size_t nb_sample) const
/* --------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. indexing file stream reference
            3. total number of sample
    Value:  1. count vector for k-mer count in query (as refArg)
            2. bool for whether the query was succeeded
            3. indexing file stream reference
    Func:   return the k-mer counts stored on disk by k-mer unique code
\* --------------------------------------------------------------------- */
{
    auto iter = kmer_serial_.find(kmer_code);
    if (iter == kmer_serial_.cend())
    {
        return false;
    }
    LoadCountFromIndex(count_vect, index_file, index_pos_.at(iter->second), nb_sample);
    return true;
}

const size_t KMerCountTab::GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set) const
/* --------------------------------------------------------------------------- *\
    Arg:    k-mer unique code set
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. bool for whether there is at least one k-mer registered
    Func:   return the average counts of given k-mer list by searching memory
\* --------------------------------------------------------------------------- */
{
    size_t nb_sample = count_tab_.at(0).size(), nb_kmer_found(0);
    std::vector<float> count_sum_vect(nb_sample, 0);
    for (auto kmer_code : kmer_code_set)
    {
        auto iter = kmer_serial_.find(kmer_code);
        if (iter != kmer_serial_.cend())
        {
            ++nb_kmer_found;
            for (unsigned int i(0); i < nb_sample; ++i)
            {
                count_sum_vect.at(i) += count_tab_.at(iter->second).at(i);
            }
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_sample; ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

const size_t KMerCountTab::GetAvgCountOnDsk(std::vector<float> &count_avg_vect,
                                            const std::set<uint64_t> kmer_code_set,
                                            std::ifstream &index_file,
                                            const size_t nb_sample) const
/* ------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code set
            2. indexing file stream reference
            3. total number of sample
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. bool for whether there is at least one k-mer registered
    Func:   return the average counts of given k-mer list by searching disk index
\* ------------------------------------------------------------------------------- */
{
    std::vector<float> count_sum_vect(nb_sample, 0);
    size_t nb_kmer_found(0);
    for (auto kmer_code : kmer_code_set)
    {
        auto iter = kmer_serial_.find(kmer_code);
        if (iter != kmer_serial_.cend())
        {
            ++nb_kmer_found;
            AddCountFromIndex(count_sum_vect, index_file, index_pos_.at(iter->second), nb_sample);
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_sample; ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}
