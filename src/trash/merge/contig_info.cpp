#include <algorithm>

#include "dna_sequence.hpp"
#include "contig_info.hpp"

ContigInfo::ContigInfo(const std::string &kmer_seq, const double rep_value, const size_t count_pos)
    : is_used_(false),
      tag_(kmer_seq),
      seq_(kmer_seq),
      nb_kmer_(1),
      rep_value_(rep_value),
      count_pos_(1, count_pos),
      head_pos_(count_pos),
      rear_pos_(count_pos)
{
}

ContigInfo::ContigInfo(const std::string &kmer_seq, const double rep_value, const std::vector<double> &count)
    : is_used_(false),
      tag_(kmer_seq),
      seq_(kmer_seq),
      nb_kmer_(1),
      rep_value_(rep_value),
      counts_(count),
      head_counts_(count),
      rear_counts_(count)
{
}

ContigInfo::ContigInfo(const ContigInfo &predtig, const ContigInfo &succtig, const size_t n_overlap, const std::string &quant_mode)
    : is_used_(false)
{
    tag_ = (predtig.GetRepValue() < succtig.GetRepValue()) ? predtig.GetTag() : succtig.GetTag();
    seq_ = predtig.GetSequence() + succtig.GetSequence().substr(n_overlap);
    nb_kmer_ = predtig.GetNbKmer() + succtig.GetNbKmer();
    rep_value_ = (predtig.GetRepValue() < succtig.GetRepValue()) ? predtig.GetRepValue() : succtig.GetRepValue();
    head_counts_ = predtig.GetHeadCounts();
    rear_counts_ = succtig.GetRearCounts();
    head_pos_ = predtig.GetHeadPos();
    rear_pos_ = succtig.GetRearPos();
    if (quant_mode == "rep")
    {
        count_pos_ = (predtig.GetRepValue() < succtig.GetRepValue()) ? predtig.GetCountPos() : succtig.GetCountPos();
        counts_ = (predtig.GetRepValue() < succtig.GetRepValue()) ? predtig.GetCounts() : succtig.GetCounts();
    }
    else if (quant_mode == "mean")
    {
        count_pos_ = predtig.GetCountPos();
        count_pos_.insert(count_pos_.end(), succtig.GetCountPos().begin(), succtig.GetCountPos().end());
        size_t nb_sample = predtig.GetCounts().size();
        for (size_t i(0); i < nb_sample; ++i)
        {
            double sum_count = predtig.GetCounts().at(i) * predtig.GetNbKmer() + succtig.GetCounts().at(i) * succtig.GetNbKmer();
            counts_.push_back(sum_count / nb_kmer_);
        }
    }
}

const bool ContigInfo::IsUsed() const
{
    return is_used_;
}

void ContigInfo::SetUsed()
{
    // std::cerr << "Set used: " << seq_ << std::endl;
    is_used_ = true;
}

const std::string &ContigInfo::GetTag() const
{
    return tag_;
}

const std::string &ContigInfo::GetSequence() const
{
    return seq_;
}

void ContigInfo::ReverseComplement()
{
    ToReverseComplement(seq_);
    head_counts_.swap(rear_counts_);
    std::swap(head_pos_, rear_pos_);
}

const size_t ContigInfo::GetNbKmer() const
{
    return nb_kmer_;
}

const double ContigInfo::GetRepValue() const
{
    return rep_value_;
}

const std::vector<size_t> &ContigInfo::GetCountPos() const
{
    return count_pos_;
}

const std::vector<double> &ContigInfo::GetCounts() const
{
    return counts_;
}

const std::vector<double> &ContigInfo::GetHeadCounts() const
{
    return head_counts_;
}

const std::vector<double> &ContigInfo::GetRearCounts() const
{
    return rear_counts_;
}

const size_t ContigInfo::GetHeadPos() const
{
    return head_pos_;
}

const size_t ContigInfo::GetRearPos() const
{
    return rear_pos_;
}
