#include "contig_elem.hpp"

ContigElem::ContigElem(const std::string &seq, const float rep_value, const uint64_t initial_kmer_code)
    : SeqElem(seq), rep_value_(rep_value)
{
    member_kmer_set_.insert(initial_kmer_code);
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

const bool ContigElem::LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap)
{
    auto left_seq = left_contig_elem.GetSeq();
    if (left_seq.size() < n_overlap || seq_.size() < n_overlap)
    {
        return false;
    }
    seq_ = left_seq + seq_.substr(n_overlap);
    member_kmer_set_.insert(left_contig_elem.GetMemberKMerSet().begin(), left_contig_elem.GetMemberKMerSet().end());
    return true;
}

const bool ContigElem::RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap)
{
    auto right_seq = right_contig_elem.GetSeq();
    if (seq_.size() < n_overlap || right_seq.size() < n_overlap)
    {
        return false;
    }
    seq_ = seq_ + right_seq.substr(n_overlap);
    member_kmer_set_.insert(right_contig_elem.GetMemberKMerSet().begin(), right_contig_elem.GetMemberKMerSet().end());
    return true;
}
