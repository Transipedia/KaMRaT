#include <sstream>
#include <iostream>

#include "sample_info.hpp"

SampleInfo::SampleInfo()
    : nb_condition_(0)
{
}

void SampleInfo::AddSampleInfo(const std::string &sample_info_line)
{
    std::istringstream conv(sample_info_line);
    std::string sample, condition;
    conv >> sample >> condition;
    if (conv.fail()) // dealing with the case when sample_info contains only sample column
    {
        condition = "_condition?absent!";
        nb_condition_ = 1;
    }
    auto iter = condition_label_.find(condition);
    size_t label;
    if (iter == condition_label_.cend())
    {
        condition_label_.insert({condition, nb_condition_});
        label = nb_condition_;
        ++nb_condition_;
    }
    else
    {
        label = iter->second;
    }
    if (!sample_label_.insert({sample, label}).second)
    {
        std::cerr << "ERROR: sample-info file has duplicated sample name (" << sample << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
}

const bool SampleInfo::IsEmpty() const
{
    return sample_label_.empty();
}

const size_t SampleInfo::GetNbCondition() const
{
    return nb_condition_;
}

const size_t SampleInfo::GetNbSample() const
{
    return sample_label_.size();
}

const int SampleInfo::GetLabel(const std::string &sample_name) const
{
    auto iter = sample_label_.find(sample_name);
    if (sample_label_.empty())
    {
        return 0;
    } 
    else if (iter == sample_label_.cend())
    {
        return -2;
    }
    else 
    {
        return iter->second;
    }
}

const bool SampleInfo::IsSample(const std::string &col_name) const
{
    auto iter = sample_label_.find(col_name);
    return (sample_label_.empty() || iter != sample_label_.cend());
}

void SampleInfo::PrintSampleNames(std::ostream &out_s, const std::string &leading_str) const
{
    out_s << leading_str;
    for (const auto &elem : sample_label_)
    {
        out_s << "\t" << elem.first;
    }
    out_s << std::endl;
}