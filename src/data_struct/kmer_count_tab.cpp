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
        if (col_info.GetColumnNature(i) == 'r' || col_info.GetColumnNature(i) == 'o')
        {
            value_vect.push_back(std::stof(term));
        }
        else if (col_info.GetColumnNature(i) == 's')
        {
            count_vect.push_back(std::stof(term));
        }
    }
}

KMerCountTab::KMerCountTab(const std::string &mode, const ColumnInfo &col_info)
    : mode_(mode), col_info_(col_info)
{
    for (unsigned int i(0); i < col_info_.GetNbColumn(); ++i)
    {
        if (col_info_.GetColumnNature(i) == 'r' || col_info_.GetColumnNature(i) == 'o')
        {
            if (!value_name_dict_.insert({col_info_.GetColumnName(i), value_name_dict_.size()}).second)
            {
                throw("duplicate column name " + col_info_.GetColumnName(i));
            }
            std::cerr << col_info_.GetColumnName(i) << "\t" << value_name_dict_.size() << std::endl;
        }
    }
}

const bool KMerCountTab::AddKMerCountInMem(const uint64_t kmer_code, const std::string &line_str)
/* -------------------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. string of input k-mer table line
    Value:  bool indicating whether the insertion was succeeded
    Func:   1. check whether the value & count vector to insert is coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. associate the k-mer unique code with its inserted row number in the value table
\* -------------------------------------------------------------------------------------------- */
{
    std::vector<float> count_vect, value_vect;
    StrLine2ValueCountVects(value_vect, count_vect, line_str, col_info_);
    if (!value_tab_.empty() && value_tab_.back().size() != value_vect.size())
    {
        throw "newly inserted k-mer value vector not coherent with existing k-mer value table";
    }
    if (!count_tab_.empty() && count_tab_.back().size() != count_vect.size())
    {
        throw "newly inserted k-mer count vector not coherent with existing k-mer count table";
    }

    bool is_succeeded = kmer_serial_.insert({kmer_code, value_tab_.size()}).second;
    if (is_succeeded)
    {
        value_tab_.emplace_back(value_vect);
        count_tab_.emplace_back(count_vect);
    }
    return is_succeeded;
}

const bool KMerCountTab::AddKMerIndexOnDsk(const uint64_t kmer_code, const std::string &line_str, std::ofstream &index_file)
/* ---------------------------------------------------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. string of input k-mer table line
            3. indexing file stream reference
    Value:  bool indicating whether the insertion was succeeded
    Func:   1. check whether the value vector to insert is coherent with the table
            2. parse the line string and index the k-mer count vector in the index file and insert value vector to value table
            3. associate the k-mer unique code with its indexed position in the file
\* ---------------------------------------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    StrLine2ValueCountVects(value_vect, count_vect, line_str, col_info_);
    if (!value_tab_.empty() && value_tab_.back().size() != value_vect.size())
    {
        throw "newly inserted k-mer value vector not coherent with existing k-mer value table";
    }
    bool is_succeeded = kmer_serial_.insert({kmer_code, value_tab_.size()}).second;
    if (is_succeeded)
    {
        value_tab_.emplace_back(value_vect);
        index_pos_.emplace_back(WriteCountToIndex(index_file, value_vect));
    }
    return is_succeeded;
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

const bool KMerCountTab::GetValue(float &value, const uint64_t kmer_code, const std::string &value_name) const
/* --------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. value name string
    Value:  1. float for the value in query (as refArg)
            2. bool for whether the query was succeeded
    Func:   return the value in query according to k-mer unique code and value name
\* --------------------------------------------------------------------------------- */
{
    auto iter_value = value_name_dict_.find(value_name);
    if (iter_value == value_name_dict_.cend())
    {
        throw "unknown value name";
    }
    auto iter_kmer = kmer_serial_.find(kmer_code);
    if (iter_kmer == kmer_serial_.cend())
    {
        return false;
    }
    value = value_tab_.at(iter_kmer->second).at(iter_value->second);
    return true;
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

const bool KMerCountTab::GetCountOnDsk(std::vector<float> &count_vect, const uint64_t kmer_code, std::ifstream &index_file) const
/* --------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code
            2. indexing file stream reference
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
    LoadCountFromIndex(count_vect, index_file, index_pos_.at(iter->second), col_info_.GetNbSample());
    return true;
}

const bool KMerCountTab::GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set) const
/* --------------------------------------------------------------------------- *\
    Arg:    k-mer unique code set
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. bool for whether there is at least one k-mer registered
    Func:   return the average counts of given k-mer list by searching memory
\* --------------------------------------------------------------------------- */
{
    std::vector<float> count_sum_vect(col_info_.GetNbSample(), 0);
    size_t nb_kmer_found(0);
    for (auto kmer_code : kmer_code_set)
    {
        auto iter = kmer_serial_.find(kmer_code);
        if (iter != kmer_serial_.cend())
        {
            ++nb_kmer_found;
            for (unsigned int i(0); i < col_info_.GetNbSample(); ++i)
            {
                count_sum_vect.at(i) += count_tab_.at(iter->second).at(i);
            }
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < col_info_.GetNbSample(); ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

const bool KMerCountTab::GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set, std::ifstream &index_file) const
/* ------------------------------------------------------------------------------- *\
    Arg:    1. k-mer unique code set
            2. indexing file stream reference
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. bool for whether there is at least one k-mer registered
    Func:   return the average counts of given k-mer list by searching disk index
\* ------------------------------------------------------------------------------- */
{
    std::vector<float> count_sum_vect(col_info_.GetNbSample(), 0);
    size_t nb_kmer_found(0);
    for (auto kmer_code : kmer_code_set)
    {
        auto iter = kmer_serial_.find(kmer_code);
        if (iter != kmer_serial_.cend())
        {
            ++nb_kmer_found;
            AddCountFromIndex(count_sum_vect, index_file, index_pos_.at(iter->second), col_info_.GetNbSample());
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < col_info_.GetNbSample(); ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}
