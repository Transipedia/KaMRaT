#include "merge_knot.hpp"

MergeKnot::MergeKnot() noexcept
    : has_pred_(false), has_succ_(false), has_ambiguity_(false)
{
}

void MergeKnot::AddContig(const size_t contig_serial, const bool is_rc, const std::string &&which_to_set) noexcept
{
    if (which_to_set == "pred" && !has_pred_)
    {
        pred_serial_ = contig_serial;
        is_pred_rc_ = is_rc;
        has_pred_ = true;
    }
    else if (which_to_set == "succ" && !has_succ_)
    {
        succ_serial_ = contig_serial;
        is_succ_rc_ = is_rc;
        has_succ_ = true;
    }
    else
    {
        has_ambiguity_ = true;
    }
}

const size_t MergeKnot::GetSerial(const std::string &&which_to_get) const noexcept
{
    if (which_to_get == "pred")
    {
        return pred_serial_;
    }
    else if (which_to_get == "succ")
    {
        return succ_serial_;
    }
    else // should not happen
    {
        return 0;
    }
    
}

const bool MergeKnot::IsRC(const std::string &&which_to_get) const noexcept
{
    if (which_to_get == "pred")
    {
        return is_pred_rc_;
    }
    else if (which_to_get == "succ")
    {
        return is_succ_rc_;
    }
    else // should not happen
    {
        return 0;
    }
    
}

const bool MergeKnot::IsMergeable() const noexcept
{
    if (has_ambiguity_)
    {
        return false;
    }
    if (!has_pred_ || !has_succ_)
    {
        return false;
    }
    if (pred_serial_ == succ_serial_)
    {
        return false;
    }
    return true;
}

const bool MergeKnot::HasPred() const noexcept
{
    return has_pred_;
}

const bool MergeKnot::HasSucc() const noexcept
{
    return has_succ_;
}