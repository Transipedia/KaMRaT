#include <stdexcept>

#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &seq, const uint64_t uniq_code, const float score)
    : seq_(seq), uniq_code_(uniq_code), score_(score), final_score_(score)
{
}

const std::string SeqElem::GetSeq() const
{
    return seq_;
}

const uint64_t SeqElem::GetUniqCode() const
{
    return uniq_code_;
}

const float SeqElem::GetScore(const std::string &&mode) const
{
    if (mode == "origin")
    {
        return score_;
    }
    else if (mode == "final")
    {
        return final_score_;
    }
    else
    {
        throw std::domain_error("unknown score mode in SeqElem class");
    }
}

void SeqElem::ScaleScore(float fact, float lower_lim, float upper_lim)
{
    final_score_ *= fact;
    final_score_ = (final_score_ < lower_lim) ? lower_lim : final_score_;
    final_score_ = (final_score_ > upper_lim) ? upper_lim : final_score_;
}