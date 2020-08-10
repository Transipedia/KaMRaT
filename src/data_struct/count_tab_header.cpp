#include <sstream>

#include "count_tab_header.hpp"

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
    size_t nb_cond(0);
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

CountTabHeader::CountTabHeader()
    : nb_cond_(0), nb_count_(0), nb_value_(0), nb_str_(0)
{
}

const void CountTabHeader::MakeColumnInfo(const std::string &header_line,
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
    colserial_vect_.push_back(nb_str_++);
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
            else if (term == "tag" || term == "feature") // specially for kamratNorm and kamratReduce
            {
                colnature_vect_.push_back('f');
                colserial_vect_.push_back(nb_str_++);
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

const std::string CountTabHeader::GetColName(const size_t i_col) const
{
    return colname_vect_.at(i_col);
}

const char CountTabHeader::GetColNature(const size_t i_col) const
{
    return colnature_vect_.at(i_col);
}

const size_t CountTabHeader::GetColSerial(const size_t i_col) const
{
    return colserial_vect_.at(i_col);
}

const size_t CountTabHeader::GetSmpLabel(const size_t i_col) const
{
    return smplabel_vect_.at(i_col);
}

const std::vector<size_t> &CountTabHeader::GetSmpLabels() const
{
    return smplabel_vect_;
}

const size_t CountTabHeader::GetNbCondition() const
{
    return nb_cond_;
}

const size_t CountTabHeader::GetNbValue() const
{
    return nb_value_;
}

const size_t CountTabHeader::GetNbCount() const
{
    return nb_count_;
}

const size_t CountTabHeader::GetNbColumn() const
{
    return colname_vect_.size();
}
