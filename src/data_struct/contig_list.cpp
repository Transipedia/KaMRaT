#include "contig_list.hpp"

ContigElem::ContigElem(const std::string &tag,
                       const std::string &seq,
                       const float rep_value,
                       const uint64_t initial_kmer_code)
    : SeqElem(tag, seq), rep_value_(rep_value)
{
    member_kmer_vect_.push_back(initial_kmer_code);
}

const float ContigElem::GetRepValue() const
{
    return rep_value_;
}

const std::vector<uint64_t> ContigElem::GetMemberKMerVect() const
{
    return member_kmer_vect_;
}

const size_t ContigElem::GetNbMemberKMer() const
{
    return member_kmer_vect_.size();
}

const bool ContigElem::LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap)
{
    if (seq_.size() < n_overlap)
    {
        return false;
    }
    seq_ = left_contig_elem.GetSeq() + seq_.substr(n_overlap);
    member_kmer_vect_.insert(member_kmer_vect_.begin(), left_contig_elem.GetMemberKMerVect().begin(), left_contig_elem.GetMemberKMerVect().end());
    return true;
}

const bool ContigElem::RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap)
{
    if (seq_.size() < n_overlap)
    {
        return false;
    }
    seq_ = seq_ + right_contig_elem.GetSeq().substr(n_overlap);
    member_kmer_vect_.insert(member_kmer_vect_.end(), right_contig_elem.GetMemberKMerVect().begin(), right_contig_elem.GetMemberKMerVect().end());
    return true;
}
