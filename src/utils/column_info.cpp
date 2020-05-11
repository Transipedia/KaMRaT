#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include "utils.hpp"
#include "column_info.hpp"

/* ------------------------------------------------------------------------------------------------ *\ 
    Args:   condition-tag dictionary (to modify), a condition string, the current condition number
    Value:  number of condition registered
    Func:   if condition does not exist in cond_tag, insert it
            return corresponding tag of the condition 
\* ------------------------------------------------------------------------------------------------ */
inline unsigned int SetConditionTag(std::unordered_map<std::string, unsigned int> &cond_tag,
                                    const std::string &condition,
                                    const unsigned int nb_cond)
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
inline unsigned int LoadSampleInfo(std::unordered_map<std::string, unsigned int> &sample_tag,
                                   const std::string &sample_info_path)

{
    std::ifstream sample_info_file(sample_info_path);
    ExitIf(!sample_info_file.is_open(), "ERROR: meta data file " + sample_info_path + " was not found.");
    unsigned int nb_cond(1);
    std::unordered_map<std::string, unsigned int> cond_tag;
    std::string sample_info_line;
    while (std::getline(sample_info_file, sample_info_line))
    {
        std::istringstream conv(sample_info_line);
        std::string sample, condition;
        conv >> sample >> condition;
        unsigned int i_tag = (!conv.fail()) ? SetConditionTag(cond_tag, condition, nb_cond) : 0;
        nb_cond = (nb_cond == i_tag) ? (i_tag + 1) : nb_cond; // cond_tag begins from 0, so nb_cond := max(i_tag) + 1
        ExitIf(!sample_tag.insert({sample, i_tag}).second, "ERROR: sample-info file has duplicated sample name " + sample + ".");
    }
    sample_info_file.close();
    return nb_cond;
}

ColumnInfo::ColumnInfo()
    : nb_sample_(0), nb_condition_(0)
{
}

void ColumnInfo::MakeColumnInfo(const std::string &header_line,
                                const std::string &sample_info_path,
                                const std::string &rep_col_name)
{
    std::unordered_map<std::string, unsigned int> sample_tag;
    nb_condition_ = (!sample_info_path.empty()) ? LoadSampleInfo(sample_tag, sample_info_path) : 1;
    std::istringstream conv(header_line);
    std::string header_term;
    // first column //
    conv >> header_term;
    ExitIf(conv.fail(), "ERROR: empty header line.");
    col_name_vect_.push_back(header_term);
    col_nat_vect_.push_back(-2); // -2 for feature column
    // other columns //
    if (sample_info_path.empty())
    {
        while (conv >> header_term)
        {
            ++nb_sample_;
            col_name_vect_.push_back(header_term);
            col_nat_vect_.push_back(0);
        }
    }
    else
    {
        while (conv >> header_term)
        {
            col_name_vect_.push_back(header_term);
            auto iter = sample_tag.find(header_term);
            if (iter != sample_tag.cend())
            {
                ++nb_sample_;
                col_nat_vect_.push_back(iter->second); // non-negative values for sample
            }
            else if (!rep_col_name.empty() && header_term == rep_col_name)
            {
                col_nat_vect_.push_back(-1); // -1 for rep_col or score_col name
            }
            else
            {
                col_nat_vect_.push_back(-3); // -3 for others
            }
        }
    }
    ExitIf(nb_sample_ == 0, "ERROR: no sample column found.");
}

const char ColumnInfo::GetColumnNature(size_t i_col) const
{
    if (col_nat_vect_.at(i_col) >= 0)
    {
        return 's'; // sample column
    }
    else if (col_nat_vect_.at(i_col) == -1)
    {
        return 'r'; // representative column
    }
    else if (col_nat_vect_.at(i_col) == -2)
    {
        return 'f'; // feature column
    }
    else if (col_nat_vect_.at(i_col) == -3)
    {
        return 'o'; // other column
    }
    else
    {
        ExitIf(true, "ERROR: unknown col_nat_code " + std::to_string(col_nat_vect_.at(i_col)));
    }
}

const unsigned int ColumnInfo::GetNbCondition() const
{
    return nb_condition_;
}

const unsigned int ColumnInfo::GetNbSample() const
{
    return nb_sample_;
}

void ColumnInfo::PrintSampleNames(std::ostream &out_s, const std::string &leading_string) const
{
    out_s << leading_string;
    size_t nb_cols(col_name_vect_.size());
    ExitIf(col_nat_vect_.size() != nb_cols, "ERROR: column names and colume nature numbers are different.");
    for (size_t i_col(0); i_col < nb_cols; ++i_col)
    {
        if (GetColumnNature(i_col) == 's')
        {
            out_s << "\t" << col_name_vect_.at(i_col);
        }
    }
    out_s << std::endl;
}