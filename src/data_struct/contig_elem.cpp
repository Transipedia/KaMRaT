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

static inline const void ReverseComplementSeq(std::string &seq)
{
    reverse(seq.begin(), seq.end());
    ToComplement(seq);
}

ContigElem::ContigElem(const std::string &seq, const size_t pos, const float val)
    : seq_(seq), rep_pos_(pos), rep_val_(val), head_pos_(pos), rear_pos_(pos)
{
}

const std::string &ContigElem::GetSeq() const
{
    return seq_;
}

const size_t ContigElem::GetRepPos() const
{
    return rep_pos_;
}

const float ContigElem::GetRepVal() const
{
    return rep_val_;
}

const size_t ContigElem::GetHeadPos(const bool need_reverse) const
{
    return (need_reverse ? rear_pos_ : head_pos_);
}

const size_t ContigElem::GetRearPos(const bool need_reverse) const
{
    return (need_reverse ? head_pos_ : rear_pos_);
}

const void ContigElem::ReverseComplement()
{
    ReverseComplementSeq(seq_);
    std::swap(head_pos_, rear_pos_);
}

const void ContigElem::LeftExtend(std::unique_ptr<ContigElem> left_contig_elem, const bool need_left_rc, unsigned int n_overlap)
{
    if (need_left_rc)
    {
        left_contig_elem->ReverseComplement();
    }
    seq_ = left_contig_elem->GetSeq() + seq_.substr(n_overlap);
    head_pos_ = left_contig_elem->GetHeadPos(false); // left contig already reversed for sequence merging
    left_contig_elem.reset();
}

const void ContigElem::RightExtend(std::unique_ptr<ContigElem> right_contig_elem, const bool need_right_rc, unsigned int n_overlap)
{
    if (need_right_rc)
    {
        right_contig_elem->ReverseComplement();
    }
    seq_ = seq_ + right_contig_elem->GetSeq().substr(n_overlap);
    rear_pos_ = right_contig_elem->GetRearPos(false); // right contig already reversed for sequence merging
    right_contig_elem.reset();
}
