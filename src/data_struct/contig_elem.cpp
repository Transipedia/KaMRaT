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
    : seq_(seq), rep_pos_(pos), rep_val_(val), mem_pos_vect_(1, pos)
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
    if (mem_pos_vect_.size() == 1)
    {
        return rep_pos_; // same as mem_pos_vect_[0]
    }
    else if (need_reverse)
    {
        return mem_pos_vect_.back();
    }
    else
    {
        return mem_pos_vect_.front();
    }
}

const size_t ContigElem::GetRearPos(const bool need_reverse) const
{
    if (mem_pos_vect_.size() == 1)
    {
        return rep_pos_; // same as mem_pos_vect_[0]
    }
    else if (need_reverse)
    {
        return mem_pos_vect_.front();
    }
    else
    {
        return mem_pos_vect_.back();
    }
}

const std::vector<size_t> &ContigElem::GetMemPosVect() const
{
    return mem_pos_vect_;
}

const size_t ContigElem::GetNbMemKmer() const
{
    return mem_pos_vect_.size();
}

const void ContigElem::ReverseComplement()
{
    ReverseComplementSeq(seq_);
    if (mem_pos_vect_.size() > 1)
    {
        std::swap(mem_pos_vect_.front(), mem_pos_vect_.back());
    }
}

const void ContigElem::LeftExtend(std::unique_ptr<ContigElem> left_contig_elem, const bool need_left_rc, unsigned int n_overlap)
{
    if (need_left_rc)
    {
        left_contig_elem->ReverseComplement();
    }
    seq_ = left_contig_elem->GetSeq() + seq_.substr(n_overlap);
    const size_t new_head_pos = mem_pos_vect_.size() - 1;
    mem_pos_vect_.insert(mem_pos_vect_.end() - 1,
                         left_contig_elem->GetMemPosVect().begin(),
                         left_contig_elem->GetMemPosVect().end()); // keep the current rear unchanged
    std::swap(mem_pos_vect_[0], mem_pos_vect_[new_head_pos]);      // assign left contig's rep pos as the new head pos
    left_contig_elem.reset();
}

const void ContigElem::RightExtend(std::unique_ptr<ContigElem> right_contig_elem, const bool need_right_rc, unsigned int n_overlap)
{
    if (need_right_rc)
    {
        right_contig_elem->ReverseComplement();
    }
    seq_ = seq_ + right_contig_elem->GetSeq().substr(n_overlap);
    mem_pos_vect_.insert(mem_pos_vect_.end(),
                         right_contig_elem->GetMemPosVect().begin(),
                         right_contig_elem->GetMemPosVect().end());
    right_contig_elem.reset();
}
