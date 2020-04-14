#include <sstream>

#include "contig_sample_info.hpp"

ContigSampleInfo::ContigSampleInfo(const std::string &contig_seq)
    : contig_seq_(contig_seq)
{
}

const std::string &ContigSampleInfo::GetContigSeq() const
{
    return contig_seq_;
}

void ContigSampleInfo::AddKmerCounts(const std::string &count_line,
                                     const std::vector<std::string> &col_name_vect,
                                     const std::vector<bool> &is_col_sample)
{
    std::istringstream conv(count_line);
    std::string str_term;
    for (int i(0); conv >> str_term; ++i)
    {
        if (is_col_sample[i])
        {
            std::string sample_name(col_name_vect[i]);
            auto iter = sample_counts_.find(sample_name);
            if (iter == sample_counts_.cend())
            {
                sample_counts_.insert({sample_name, std::vector<size_t>(1, std::stod(str_term) + 0.5)}); //add 0.5 for rounding
            }
            else
            {
                iter->second.push_back(std::stod(str_term) + 0.5); //add 0.5 for rounding
            }
        }
    }
}

void ContigSampleInfo::EstimatePrint(std::ostream &out_s,
                                     double (*eval_func)(const std::vector<size_t> &counts),
                                     const std::vector<std::string> &col_name_vect,
                                     const std::vector<bool> &is_col_sample) const
{
    out_s << contig_seq_;
    for (int i(0); i < col_name_vect.size(); ++i)
    {
        if (is_col_sample[i])
        {
            std::string sample_name(col_name_vect[i]);
            out_s << "\t" << eval_func(sample_counts_.find(sample_name)->second);
        }
    }
    out_s << std::endl;
}
