#include <algorithm>

#include "contig_elem.hpp"

inline void ToComplement(std::string &seq)
{
    auto lambda_trans = [](const char c) {
        switch (c)
        {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            return 'N';
        }
    };
    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda_trans);
}

ContigElem::ContigElem(const std::string &seq, const float rep_value, const size_t init_serial)
    : SeqElem(seq),
      is_used_(false),
      rep_value_(rep_value),
      rep_serial_(init_serial),
      head_serial_(init_serial),
      rear_serial_(init_serial)
{
    kmer_serial_set_.insert(init_serial);
}

const size_t ContigElem::GetHeadSerial() const
{
    return head_serial_;
}

const size_t ContigElem::GetRearSerial() const
{
    return rear_serial_;
}

const size_t ContigElem::GetRepSerial() const
{
    return rep_serial_;
}

const float ContigElem::GetRepValue() const
{
    return rep_value_;
}

const std::set<size_t> &ContigElem::GetKMerSerialSet() const
{
    return kmer_serial_set_;
}

const size_t ContigElem::GetNbMemberKMer() const
{
    return kmer_serial_set_.size();
}

const void ContigElem::LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap)
{
    auto left_seq = left_contig_elem.GetSeq();
    seq_ = left_seq + seq_.substr(n_overlap);
    kmer_serial_set_.insert(left_contig_elem.GetKMerSerialSet().begin(), left_contig_elem.GetKMerSerialSet().end());
    head_serial_ = left_contig_elem.GetHeadSerial();
}

const void ContigElem::RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap)
{
    auto right_seq = right_contig_elem.GetSeq();
    seq_ = seq_ + right_seq.substr(n_overlap);
    kmer_serial_set_.insert(right_contig_elem.GetKMerSerialSet().begin(), right_contig_elem.GetKMerSerialSet().end());
    rear_serial_ = right_contig_elem.GetRearSerial();
}

const void ContigElem::SelfReverseComplement()
{
    reverse(seq_.begin(), seq_.end());
    ToComplement(seq_);
    std::swap(head_serial_, rear_serial_);
}
