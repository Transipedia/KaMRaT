#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include "column_info.hpp"

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

ColumnInfo::ColumnInfo()
    : nb_condition_(0), nb_sample_(0), nb_value_(0)
{
}

void ColumnInfo::MakeColumnInfo(const std::string &header_line,
                                const std::string &sample_info_path)
{
    if (header_line.empty()) // quit if header line is empty
    {
        throw std::domain_error("cannot make ColumnInfo object with an empty header line");
    }

    //----- Create Sample-Tag Dictionary -----//
    std::map<std::string, size_t> sample_tag;
    nb_condition_ = (!sample_info_path.empty()) ? LoadSampleInfo(sample_tag, sample_info_path) : 1;

    //----- Make ColumnInfo Object -----//
    std::istringstream conv(header_line);
    std::string term;
    conv >> term; // the first term which is supposed to be always the tag/feature column
    colname_vect_.push_back(term);
    colname2label_.insert({term, -2});
    colname2serial_.insert({term, 0});

    if (sample_info_path.empty()) // if no sample info file provided => all next columns are samples
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            colname2label_.insert({term, 0});
            colname2serial_.insert({term, nb_sample_++});
        }
    }
    else // if a sample info file is provided => next columns may contain non-sample values
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            auto iter = sample_tag.find(term);
            if (iter != sample_tag.cend())
            {
                colname2label_.insert({term, iter->second}); // non-negative values for sample
                colname2serial_.insert({term, nb_sample_++});
            }
            else
            {
                colname2label_.insert({term, -1}); // -1 for non-count value
                colname2serial_.insert({term, nb_value_++});
            }
        }
    }
    if (nb_sample_ == 0)
    {
        throw std::domain_error("no sample column found");
    }
}

const std::string ColumnInfo::GetColumnName(const size_t i_col) const
{
    return colname_vect_.at(i_col);
}

const int ColumnInfo::GetColumnLabel(const size_t i_col) const
{
    return colname2label_.find(colname_vect_.at(i_col))->second;
}

const char ColumnInfo::GetColumnNature(const size_t i_col) const
{
    if (colname2label_.find(colname_vect_.at(i_col))->second >= 0)
    {
        return 's'; // sample-count column
    }
    else if (colname2label_.find(colname_vect_.at(i_col))->second == -1)
    {
        return 'v'; // value column
    }
    else if (colname2label_.find(colname_vect_.at(i_col))->second == -2)
    {
        return 'f'; // feature column
    }
    else
    {
        throw std::domain_error("unknown column nature code " + std::to_string(colname2label_.find(colname_vect_.at(i_col))->second));
    }
}

const size_t ColumnInfo::GetColumnSerial(const std::string &colname) const
{
    auto iter = colname2serial_.find(colname);
    if (iter == colname2serial_.cend())
    {
        throw std::domain_error("non-registered column name " + colname);
    }
    return iter->second;
}

const size_t ColumnInfo::GetNbCondition() const
{
    return nb_condition_;
}

const size_t ColumnInfo::GetNbSample() const
{
    return nb_sample_;
}

const size_t ColumnInfo::GetNbColumn() const
{
    return colname_vect_.size();
}