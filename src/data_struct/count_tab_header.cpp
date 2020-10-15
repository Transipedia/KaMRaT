#include <sstream>

#include "count_tab_header.hpp"

inline size_t StrLine2ValueCountVects(std::vector<float> &value_vect,
                                      std::vector<float> &count_vect,
                                      std::vector<std::string> &str_vect,
                                      const std::string &line_str,
                                      const std::vector<int> &colnature_vect)
{
    std::istringstream conv(line_str);
    std::string term;
    size_t score_pos(0); // ASSUMPTION: score value never appears at the first column !!!
    conv >> term;        // skip the first column which is k-mer/tag/contig/sequence
    for (unsigned int i(1); conv >> term && i < colnature_vect.size(); ++i)
    {
        if (colnature_vect[i] == -2)
        {
            value_vect.push_back(std::stof(term));
        }
        else if (colnature_vect[i] == -1)
        {
            score_pos = i;
            value_vect.push_back(std::stof(term));
        }
        else if (colnature_vect[i] >= 0)
        {
            count_vect.push_back(std::stof(term));
        }
        else if (colnature_vect[i] == -3)
        {
            str_vect.push_back(term);
        }
        else
        {
            throw std::invalid_argument("unknown column nature code: " + std::to_string(colnature_vect[i]));
        }
    }
    if (value_vect.size() + count_vect.size() + str_vect.size() != colnature_vect.size() - 1) // the first column is excluded for saving memory
    {
        throw std::domain_error("parsing string line to fields failed: sum of field lengths not equal to header length");
    }
    return score_pos;
}

CountTabHeader::CountTabHeader()
    : nb_col_(0), nb_cond_(0), nb_value_(0), nb_count_(0), nb_str_(0)
{
}

const void CountTabHeader::MakeSmpCond(const std::string &sample_info_path,
                                       const std::set<std::string> &preserved_cond_tags)
{
    if (!sample_info_path.empty())
    {
        std::ifstream sample_info_file(sample_info_path);
        if (!sample_info_file.is_open())
        {
            throw std::domain_error("meta data file " + sample_info_path + " was not found");
        }
        std::string sample_info_line;
        while (std::getline(sample_info_file, sample_info_line))
        {
            std::istringstream conv(sample_info_line);
            std::string sample, condition;
            conv >> sample >> condition;
            if (conv.fail())
            {
                condition = ""; // empty string indicating the unique condition
            }
            if (preserved_cond_tags.find(condition) != preserved_cond_tags.cend())
            {
                throw std::domain_error("condition name " + condition + " is preserved by KaMRaT, please change another name for this condition");
            }
            auto ins = cond2lab_dict_.insert({condition, nb_cond_});
            size_t i_tag = ins.first->second;
            if (!smp2lab_dict_.insert({sample, i_tag}).second)
            {
                throw std::domain_error("sample-info file has duplicated sample name " + sample);
            }
            if (nb_cond_ == i_tag)
            {
                ++nb_cond_;
            }
        }
        sample_info_file.close();
    }
}

const void CountTabHeader::MakeSmpCond(const std::string &sample_info_path)
{
    const std::set<std::string> preserved_cond_tags;
    MakeSmpCond(sample_info_path, preserved_cond_tags);
}

const void CountTabHeader::MakeColumnInfo(const std::string &header_line,
                                          const std::string &score_colname)
{
    if (header_line.empty()) // quit if header line is empty
    {
        throw std::domain_error("cannot parse column information with an empty header line");
    }
    std::istringstream conv(header_line);
    std::string term;
    // first term is supposed to be feature //
    conv >> term;
    colname_vect_.push_back(term);
    colnature_vect_.push_back(-3);
    colserial_vect_.push_back(nb_str_++);
    ++nb_col_;
    // following columns //
    if (cond2lab_dict_.empty()) // sample info file NOT provided => all next columns are samples
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            colnature_vect_.push_back(0);
            colserial_vect_.push_back(nb_count_++);
            ++nb_col_;
        }
    }
    else // sample info file provided => next columns may contain non-sample values
    {
        while (conv >> term)
        {
            colname_vect_.push_back(term);
            auto iter = smp2lab_dict_.find(term);
            if (iter != smp2lab_dict_.cend())
            {
                colnature_vect_.push_back(iter->second);
                colserial_vect_.push_back(nb_count_++);
            }
            else if (term == "tag" || term == "feature" || term == "contig")
            {
                colnature_vect_.push_back(-3); // -3 for strings
                colserial_vect_.push_back(nb_str_++);
            }
            else
            {
                size_t nat_code = (term == score_colname ? -1 : -2); // -1 for rep-value or score, -2 for ordinary value
                colnature_vect_.push_back(nat_code);
                colserial_vect_.push_back(nb_value_++);
            }
            ++nb_col_;
        }
    }
    if (nb_count_ == 0)
    {
        throw std::domain_error("no sample column found");
    }
    if (nb_col_ != colname_vect_.size() || nb_col_ != colserial_vect_.size() || nb_col_ != colnature_vect_.size())
    {
        throw std::domain_error("colname, colserial, colnature sizes not coherent");
    }
}

const size_t CountTabHeader::GetCondLabel(const std::string &cond_name) const
{
    auto iter = cond2lab_dict_.find(cond_name);
    if (iter == cond2lab_dict_.cend())
    {
        throw std::domain_error("condition name not exists: " + cond_name);
    }
    return iter->second;
}

const std::string CountTabHeader::GetColName(const size_t i_col) const
{
    return colname_vect_.at(i_col);
}

const int CountTabHeader::GetColNature(const size_t i_col) const
{
    return colnature_vect_.at(i_col);
}

const size_t CountTabHeader::GetColSerial(const size_t i_col) const
{
    return colserial_vect_.at(i_col);
}

const void CountTabHeader::GetSmpLabels(std::vector<size_t> &smp_labels)
{
    if (!smp_labels.empty())
    {
        throw std::domain_error("sample label vector not empty before adding sample labels");
    }
    for (size_t i(0); i < nb_col_; ++i)
    {
        int colnat = colnature_vect_[i];
        if (colnat >= 0)
        {
            smp_labels.push_back(colnat);
        }
    }
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
    return nb_col_;
}

const bool CountTabHeader::IsSample(size_t i_col) const
{
    return (colnature_vect_[i_col] >= 0);
}

const bool CountTabHeader::CheckColName(size_t i_col, const std::string &target_colname) const
{
    return (colname_vect_[i_col] == target_colname);
}

const bool CountTabHeader::CheckColCondition(size_t i_col, const std::string &target_condition) const
{
    auto iter = cond2lab_dict_.find(target_condition);
    if (iter == cond2lab_dict_.cend())
    {
        throw std::domain_error("unknown target condition: " + target_condition);
    }
    if (colnature_vect_[i_col] < 0)
    {
        throw std::domain_error("checking a non-sample column is not possible");
    }
    return (iter->second == static_cast<size_t>(colnature_vect_[i_col]));
}

const size_t CountTabHeader::ParseLineStr(std::vector<float> &value_vect,
                                          std::vector<float> &count_vect,
                                          std::vector<std::string> &str_vect,
                                          const std::string &line_str) const
{
    size_t score_pos = StrLine2ValueCountVects(value_vect, count_vect, str_vect, line_str, colnature_vect_);
    if (nb_value_ != value_vect.size())
    {
        throw std::domain_error("newly inserted value vector not coherent with existing value table");
    }
    if (nb_count_ != count_vect.size())
    {
        throw std::domain_error("newly inserted count vector not coherent with existing count table");
    }
    if (nb_str_ != str_vect.size() + 1) // the first column is excluded for saving memory
    {
        throw std::domain_error("newly inserted string vector not coherent with existing string table");
    }
    return score_pos;
}
