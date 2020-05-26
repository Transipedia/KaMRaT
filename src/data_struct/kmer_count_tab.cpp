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

/* ------------------------------------------------------------------------------------------------ *\ 
    Args:   condition-tag dictionary (to modify), a condition string, the current condition number
    Value:  number of condition registered
    Func:   if condition does not exist in cond_tag, insert it
            return corresponding tag of the condition 
\* ------------------------------------------------------------------------------------------------ */
inline size_t SetConditionTag(std::map<std::string, size_t> &cond_tag,
                              const std::string &condition,
                              const size_t nb_cond)
{
    auto iter = cond_tag.find(condition);
    if (iter == cond_tag.cend())
    {
        cond_tag.insert({condition, nb_cond});
        return nb_cond;
    }
    else
    {
        return iter->second;
    }
}

/* ------------------------------------------------------------------ *\
    Args:   sample-tag dictionary (to modify), sample-info file path 
    Value:  total number of registed conditions
\* ------------------------------------------------------------------ */
inline size_t LoadSampleInfo(std::map<std::string, size_t> &sample_tag,
                             const std::string &sample_info_path)

{
    std::ifstream sample_info_file(sample_info_path);
    if (!sample_info_file.is_open())
    {
        throw std::domain_error("meta data file " + sample_info_path + " was not found");
    }
    size_t nb_cond(1);
    std::map<std::string, size_t> cond_tag;
    std::string sample_info_line;
    while (std::getline(sample_info_file, sample_info_line))
    {
        std::istringstream conv(sample_info_line);
        std::string sample, condition;
        conv >> sample >> condition;
        size_t i_tag = (!conv.fail()) ? SetConditionTag(cond_tag, condition, nb_cond) : 0;
        nb_cond = (nb_cond == i_tag) ? (i_tag + 1) : nb_cond; // cond_tag begins from 0, so nb_cond := max(i_tag) + 1
        if (!sample_tag.insert({sample, i_tag}).second)
        {
            throw std::domain_error("sample-info file has duplicated sample name " + sample);
        }
    }
    sample_info_file.close();
    return nb_cond;
}

inline size_t StrLine2ValueCountVects(std::vector<float> &value_vect,
                                      std::vector<float> &count_vect,
                                      const std::string &line_str,
                                      const std::vector<char> &colnature_vect)
{
    std::istringstream conv(line_str);
    std::string term;
    size_t score_pos(0); // ASSUMPTION: score value never appears at the first column !!!
    for (unsigned int i(0); conv >> term; ++i)
    {
        if (colnature_vect.at(i) == 'v')
        {
            value_vect.push_back(std::stof(term));
        }
        else if (colnature_vect.at(i) == '+')
        {
            score_pos = i;
            value_vect.push_back(std::stof(term));
        }
        else if (colnature_vect.at(i) == 's')
        {
            count_vect.push_back(std::stof(term));
        }
    }
    return score_pos;
}

KMerCountTab::KMerCountTab(const std::string &mode)
    : mode_(mode), nb_cond_(0), nb_count_(0), nb_value_(0)
{
}

const std::string KMerCountTab::GetMode() const
{
    return mode_;
}

const void KMerCountTab::MakeColumnInfo(const std::string &header_line,
                                        const std::string &sample_info_path,
                                        const std::string &score_colname)
{
    if (header_line.empty()) // quit if header line is empty
    {
        throw std::domain_error("cannot make ColumnInfo object with an empty header line");
    }

    //----- Create Sample-Tag Dictionary -----//
    std::map<std::string, size_t> sample_tag;
    nb_cond_ = (!sample_info_path.empty()) ? LoadSampleInfo(sample_tag, sample_info_path) : 1;

    //----- Make ColumnInfo Object -----//
    std::istringstream conv(header_line);
    std::string term;
    // first term is supposed to be feature //
    conv >> term;
    colname_vect_.push_back(term);
    colnature_vect_.push_back('f');
    colserial_vect_.push_back(0);
    // following columns //
    if (sample_info_path.empty()) // sample info file NOT provided => all next columns are samples
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            colnature_vect_.push_back('s');
            smplabel_vect_.push_back(0);
            colserial_vect_.push_back(nb_count_++);
        }
    }
    else // sample info file provided => next columns may contain non-sample values
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            auto iter = sample_tag.find(term);
            if (iter != sample_tag.cend())
            {
                colnature_vect_.push_back('s');
                smplabel_vect_.push_back(iter->second); // non-negative values for sample
                colserial_vect_.push_back(nb_count_++);
            }
            else
            {
                char nat = (term == score_colname ? '+' : 'v');
                colnature_vect_.push_back(nat);
                colserial_vect_.push_back(nb_value_++);
            }
        }
    }
    if (nb_count_ == 0)
    {
        throw std::domain_error("no sample column found");
    }
}

const std::string KMerCountTab::GetColName(const size_t i_col) const
{
    return colname_vect_.at(i_col);
}

const char KMerCountTab::GetColNature(const size_t i_col) const
{
    return colnature_vect_.at(i_col);
}

const size_t KMerCountTab::GetColSerial(const size_t i_col) const
{
    return colserial_vect_.at(i_col);
}

const int KMerCountTab::GetSmpLabel(const size_t i_col) const
{
    return smplabel_vect_.at(i_col);
}

const size_t KMerCountTab::GetNbCondition() const
{
    return nb_cond_;
}

const size_t KMerCountTab::GetNbValue() const
{
    return nb_value_;
}

const size_t KMerCountTab::GetNbCount() const
{
    return nb_count_;
}

const size_t KMerCountTab::GetNbColumn() const
{
    return colname_vect_.size();
}

const float KMerCountTab::AddKMerCountInMem(const std::string &line_str)
/* ------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. return k-mer's corresponded score value
\* ------------------------------------------------------------------------------------------- */
{
    std::vector<float> count_vect, value_vect;
    size_t score_pos = StrLine2ValueCountVects(value_vect, count_vect, line_str, colnature_vect_);
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    if (nb_value_ > 0)
    {
        value_tab_.emplace_back(value_vect);
    }
    count_tab_.emplace_back(count_vect);
    return (score_pos == 0 ? (count_tab_.size() - 1) : value_vect.at(score_pos));
}

const float KMerCountTab::AddKMerIndexOnDsk(const std::string &line_str, std::ofstream &index_file)
/* -------------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
            3. indexing file stream reference
    Value:  inserted k-mer score (a certain column's value or the input serial number)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. index the k-mer count vector in the index file and insert value vector to value table
            4. return k-mer's corresponded score value
\* -------------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect, count_vect;
    size_t score_pos = StrLine2ValueCountVects(value_vect, count_vect, line_str, colnature_vect_);
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted k-mer value vector not coherent with existing k-mer value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted k-mer count vector not coherent with existing k-mer count table");
    }
    if (nb_value_ > 0)
    {
        value_tab_.emplace_back(value_vect);
    }
    index_pos_.emplace_back(WriteCountToIndex(index_file, count_vect));
    return (score_pos == 0 ? (count_tab_.size() - 1) : value_vect.at(score_pos));
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
