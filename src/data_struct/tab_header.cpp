#include <iostream>

#include "tab_header.hpp"

TabHeader::TabHeader()
    : nb_count_(0), rep_colpos_(0), nb_condi_(1)
{
}

TabHeader::TabHeader(const std::string &smp_info_path)
    : nb_count_(0), rep_colpos_(0)
{
    std::ifstream smp_info_file(smp_info_path);
    if (!smp_info_file.is_open())
    {
        throw std::domain_error("error open sample-info file: " + smp_info_path);
    }
    std::string sample_info_line, sample, condition;
    std::istringstream conv;
    nb_condi_ = 0;
    while (std::getline(smp_info_file, sample_info_line))
    {
        conv.str(sample_info_line);
        conv >> sample >> condition;
        if (conv.fail())
        {
            condition = ""; // empty string indicating the unique condition
        }
        const auto &ins_pair = condi2lab_.insert({condition, nb_condi_});
        if (ins_pair.second) // if a new condition is added to the dictionary
        {
            ++nb_condi_;
        }
        if (!smp2lab_.insert({sample, ins_pair.first->second}).second) // associate the sample with its condition label
        {
            throw std::domain_error("sample-info file has duplicated sample name: " + sample);
        }
        conv.clear();
    }
    smp_info_file.close();
    if (nb_condi_ != condi2lab_.size()) // should not happen, for debug
    {
        throw std::domain_error("number of condition not consistent");
    }
}

const void TabHeader::MakeColumnInfo(std::istringstream &line_conv, const std::string &rep_colname)
{
    std::string term;
    // the first column is supposed to be feature string //
    line_conv >> term;
    colname_vect_.emplace_back(term);
    is_count_.emplace_back(false);
    // following columns //
    if (smp2lab_.empty())
    {
        while (line_conv >> term)
        {
            colname_vect_.emplace_back(term);
            is_count_.emplace_back(true);
            sample_labels_.emplace_back(0);
            ++nb_count_;
        }
    }
    else
    {
        while (line_conv >> term)
        {
            colname_vect_.emplace_back(term);
            const auto &iter = smp2lab_.find(term);
            if (iter != smp2lab_.cend())
            {
                is_count_.emplace_back(true);
                sample_labels_.emplace_back(iter->second);
                ++nb_count_;
            }
            else
            {
                if (rep_colname == term)
                {
                    rep_colpos_ = colname_vect_.size() - 1; // minus 1 for the current column position
                }
                is_count_.emplace_back(false);
            }
        }
    }
    if (is_count_.size() != colname_vect_.size() || sample_labels_.size() != nb_count_) // should not happen, for debug
    {
        throw std::domain_error("colname and column condition sizes not coherent");
    }
}

const std::string &TabHeader::MakeOutputHeaderStr(std::string &header) const
{
    header.clear();
    header = colname_vect_[0];
    for (size_t i(1); i < colname_vect_.size(); ++i)
    {
        if (!IsCount(i))
        {
            header += ("\t" + colname_vect_[i]);
        }
    }
    for (size_t i(1); i < colname_vect_.size(); ++i)
    {
        if (IsCount(i))
        {
            header += ("\t" + colname_vect_[i]);
        }
    }
    return header;
}

const size_t TabHeader::GetNbCondition() const
{
    return nb_condi_;
}

const size_t TabHeader::GetConditionLabel(const std::string &condi) const
{
    auto iter = condi2lab_.find(condi);
    return (iter == condi2lab_.cend() ? iter->second : '\0');
}

const size_t TabHeader::GetNbCol() const
{
    return colname_vect_.size();
}

const size_t TabHeader::GetNbCount() const
{
    return nb_count_;
}

const size_t TabHeader::GetRepColPos() const
{
    return rep_colpos_;
}

const std::string &TabHeader::GetColNameAt(const size_t i) const
{
    return colname_vect_[i];
}

const bool TabHeader::IsCount(const size_t i) const
{
    return is_count_[i];
}

const double TabHeader::ParseRowStr(std::vector<float> &count_vect, std::string &non_count_str, std::istringstream &line_conv) const
{
    count_vect.reserve(nb_count_);
    double rep_val(0);
    static std::string term;
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    non_count_str += term;
    for (size_t i(1); (i < colname_vect_.size()) && (line_conv >> term); ++i)
    {
        if (IsCount(i))
        {
            count_vect.emplace_back(std::stof(term));
        }
        else
        {
            non_count_str += ("\t" + term);
            if (i == rep_colpos_)
            {
                rep_val = std::stod(term);
            }
        }
    }
    if (count_vect.size() != nb_count_) // should not happen, for debug
    {
        throw std::domain_error("parsing string line to fields failed");
    }
    return rep_val;
}

const std::vector<size_t> &TabHeader::GetSmpLabels() const
{
    return sample_labels_;
}

const void TabHeader::PrintSmp2Lab() const
{
    for (const auto &p : smp2lab_)
    {
        std::cerr << p.first << "\t" << p.second << std::endl;
    }
}
