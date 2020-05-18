#include "contig_elem.hpp"

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

