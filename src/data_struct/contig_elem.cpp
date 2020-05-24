#include "contig_elem.hpp"

ContigElem::ContigElem(const std::string &seq, const float rep_value, const uint64_t initial_kmer_code)
    : SeqElem(seq), is_used_(false), rep_value_(rep_value)
{
    member_kmer_set_.insert(initial_kmer_code);
}

const bool ContigElem::IsUsed() const
{
    return is_used_;
}

const void ContigElem::SetUsed()
{
    is_used_ = true;
}

const float ContigElem::GetRepValue() const
{
    return rep_value_;
}

const std::set<uint64_t> ContigElem::GetMemberKMerSet() const
{
    return member_kmer_set_;
}

const size_t ContigElem::GetNbMemberKMer() const
{
    return member_kmer_set_.size();
}

const bool ContigElem::LeftMerge(ContigElem &left_contig_elem, unsigned int n_overlap)
{
    auto left_seq = left_contig_elem.GetSeq();
    if (left_seq.size() < n_overlap || seq_.size() < n_overlap)
    {
        return false;
    }
    seq_ = left_seq + seq_.substr(n_overlap);
    member_kmer_set_.insert(left_contig_elem.GetMemberKMerSet().begin(), left_contig_elem.GetMemberKMerSet().end());
    left_contig_elem.SetUsed();
    return true;
}

const bool ContigElem::RightMerge(ContigElem &right_contig_elem, unsigned int n_overlap)
{
    auto right_seq = right_contig_elem.GetSeq();
    if (seq_.size() < n_overlap || right_seq.size() < n_overlap)
    {
        return false;
    }
    seq_ = seq_ + right_seq.substr(n_overlap);
    member_kmer_set_.insert(right_contig_elem.GetMemberKMerSet().begin(), right_contig_elem.GetMemberKMerSet().end());
    right_contig_elem.SetUsed();
    return true;
}
