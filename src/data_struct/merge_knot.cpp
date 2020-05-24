#include "merge_knot.hpp"

MergeKnot::MergeKnot()
    : has_pred_(false), has_succ_(false), has_ambiguity_(false)
{
}

void MergeKnot::AddContig(const uint64_t contig_code, const bool is_rc, const std::string &which_to_set)
{
    if (which_to_set == "pred" && !has_pred_)
    {
        pred_code_ = contig_code;
        is_pred_rc_ = is_rc;
        has_pred_ = true;
    }
    else if (which_to_set == "pred")
    {
        has_ambiguity_ = true;
    }
    else if (which_to_set == "succ" && !has_succ_)
    {
        succ_code_ = contig_code;
        is_succ_rc_ = is_rc;
        has_succ_ = true;
    }
    else if (which_to_set == "succ")
    {
        has_ambiguity_ = true;
    }
    else
    {
        throw "unknown argument which_to_set " + which_to_set;
    }
}

const uint64_t MergeKnot::GetPredCode() const
{
    return pred_code_;
}

const uint64_t MergeKnot::GetSuccCode() const
{
    return succ_code_;
}

const bool MergeKnot::IsPredRC() const
{
    return is_pred_rc_;
}

const bool MergeKnot::IsSuccRC() const
{
    return is_succ_rc_;
}

const bool MergeKnot::IsMergeable() const
{
    if (has_ambiguity_)
    {
        return false;
    }
    if (!has_pred_ || !has_succ_)
    {
        return false;
    }
    if (pred_code_ == succ_code_)
    {
        return false;
    }
    return true;
}