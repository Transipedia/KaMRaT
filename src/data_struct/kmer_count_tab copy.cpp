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

inline size_t StrLine2ValueCountVects(std::vector<float> &value_vect,
                                      std::vector<float> &count_vect,
                                      const std::string &line_str,
                                      const ColumnInfo &col_info,
                                      const std::string &score_colname)
{
    std::istringstream conv(line_str);
    std::string term;
    size_t score_pos(0); // ASSUMPTION: score value never appears at the first column !!!
    for (unsigned int i(0); conv >> term; ++i)
    {
        if (col_info.GetColNature(i) == 'v')
        {
            value_vect.push_back(std::stof(term));
        }
        else if (col_info.GetColNature(i) == '+')
        {
            score_pos = i;
            value_vect.push_back(std::stof(term));
        }
        else if (col_info.GetColNature(i) == 's')
        {
            count_vect.push_back(std::stof(term));
        }
    }
    return score_pos;
}

KMerCountTab::KMerCountTab(const std::string &mode)
    : mode_(mode), nb_count_(0), nb_value_(0)
{
}

const std::string KMerCountTab::GetMode() const
{
    return mode_;
}

const float KMerCountTab::AddKMerCountInMem(const std::string &line_str,
                                            const ColumnInfo &column_info,
                                            const std::string &score_colname)
/* ------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. ColumnInfo object
            3. column name for score
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. return k-mer's corresponded score value
\* ------------------------------------------------------------------------------------------- */
{
    std::vector<float> count_vect, value_vect;
    size_t score_pos = StrLine2ValueCountVects(value_vect, count_vect, line_str, column_info, score_colname);
    if (count_tab_.empty())
    {
        nb_count_ = column_info.GetNbCount();
        nb_value_ = column_info.GetNbValue();
    }
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    value_tab_.emplace_back(value_vect);
    count_tab_.emplace_back(count_vect);
    return (score_pos == 0 ? value_tab_.size() : value_vect.at(score_pos));
}

const float KMerCountTab::AddKMerIndexOnDsk(const std::string &line_str,
                                            const ColumnInfo &column_info,
                                            const std::string &score_colname,
                                            std::ofstream &index_file)
/* -------------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. ColumnInfo object
            3. column name for score
            4. indexing file stream reference
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. index the k-mer count vector in the index file and insert value vector to value table
            4. return k-mer's corresponded score value
\* -------------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    size_t score_pos = StrLine2ValueCountVects(value_vect, count_vect, line_str, column_info, score_colname);
    if (count_tab_.empty())
    {
        nb_count_ = column_info.GetNbCount();
        nb_value_ = column_info.GetNbValue();
    }
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    value_tab_.emplace_back(value_vect);
    index_pos_.emplace_back(WriteCountToIndex(index_file, count_vect));
    return (score_pos == 0 ? value_tab_.size() : value_vect.at(score_pos));
}

const float KMerCountTab::GetValue(const size_t kmer_serial, const size_t valcol_serial) const
/* ------------------------------------------------------------------------------ *\
    Arg:    1. k-mer serial
            2. value column serial number
    Value:  float for the value in query
    Func:   return the value in query according to k-mer and value column serial
\* ------------------------------------------------------------------------------ */
{
    return value_tab_.at(kmer_serial).at(valcol_serial);
}

const void KMerCountTab::GetCountInMem(std::vector<float> &count_vect, const size_t kmer_serial) const
/* ------------------------------------------------------------------ *\
    Arg:    k-mer serial
    Value:  count vector for k-mer count in query (as refArg)
    Func:   return the k-mer counts stored in memory by k-mer serial
\* ------------------------------------------------------------------ */
{
    count_vect = count_tab_.at(kmer_serial);
}

const void KMerCountTab::GetCountOnDsk(std::vector<float> &count_vect, const size_t kmer_serial, std::ifstream &index_file) const
/* --------------------------------------------------------------------- *\
    Arg:    1. k-mer serial
            2. indexing file stream reference
    Value:  count vector for k-mer count in query (as refArg)
    Func:   return the k-mer counts stored on disk by k-mer unique code
\* --------------------------------------------------------------------- */
{
    LoadCountFromIndex(count_vect, index_file, index_pos_.at(kmer_serial), nb_count_);
}

const size_t KMerCountTab::GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set) const
/* ------------------------------------------------------------------------------------------------------ *\
    Arg:    k-mer serial set
    Value:  1. average count vector for k-mer count in query (as refArg)
            2. size_t for number of k-mer found in the table
    Func:   return the number of found k-mers and average counts of given k-mer list by searching memory
\* ------------------------------------------------------------------------------------------------------ */
{
    size_t nb_kmer_found(0);
    std::vector<float> count_sum_vect(nb_count_, 0);
    for (auto kmer_serial : kmer_serial_set)
    {
        ++nb_kmer_found;
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_sum_vect.at(i) += count_tab_.at(kmer_serial).at(i);
        }
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

const size_t KMerCountTab::GetAvgCountOnDsk(std::vector<float> &count_avg_vect,
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
    std::vector<float> count_sum_vect(nb_count_, 0);
    size_t nb_kmer_found(0);
    for (auto kmer_serial : kmer_serial_set)
    {
        ++nb_kmer_found;
        AddCountFromIndex(count_sum_vect, index_file, index_pos_.at(kmer_serial), nb_count_);
    }
    if (nb_kmer_found > 0)
    {
        for (unsigned int i(0); i < nb_count_; ++i)
        {
            count_sum_vect.at(i) /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}
