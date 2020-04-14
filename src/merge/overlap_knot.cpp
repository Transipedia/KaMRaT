#include "overlap_knot.hpp"

OverlapKnot::OverlapKnot()
    : predtig_serial_(0),
      succtig_serial_(0),
      predtig_rc_(false),
      succtig_rc_(false),
      has_ambiguity_(false)
{
}

void OverlapKnot::SetPredtig(const size_t contig_serial, const bool is_rc)
{
    if (predtig_serial_ == 0) // if the predtig has not been set
    {
        predtig_serial_ = contig_serial;
        predtig_rc_ = is_rc;
        // std::cerr << predtig_serial_ << " == " << succtig_serial_ << std::endl;
    }
    else
    {
        // std::cerr << "predtig ambiguity:\t" << predtig_serial_ << "\t" << contig_serial << " == " << succtig_serial_ << std::endl;
        has_ambiguity_ = true;
    }
}

void OverlapKnot::SetSucctig(const size_t contig_serial, const bool is_rc)
{
    if (succtig_serial_ == 0) // if the succtig has not been set
    {
        succtig_serial_ = contig_serial;
        succtig_rc_ = is_rc;
        // std::cerr << predtig_serial_ << " == " << succtig_serial_ << std::endl;
    }
    else
    {
        // std::cerr << "succtig ambiguity:\t" << predtig_serial_ << " == " << succtig_serial_ << "\t" << contig_serial << std::endl;
        has_ambiguity_ = true;
    }
}

const size_t OverlapKnot::GetPredtigSerial() const
{
    return predtig_serial_;
}

const size_t OverlapKnot::GetSucctigSerial() const
{
    return succtig_serial_;
}

const bool OverlapKnot::IsPredtigRC() const
{
    return predtig_rc_;
}

const bool OverlapKnot::IsSucctigRC() const
{
    return succtig_rc_;
}

const bool OverlapKnot::IsSingleKnot() const
{
    if (has_ambiguity_)
    {
        return false;
    }
    if (predtig_serial_ == 0 || succtig_serial_ == 0)
    {
        return false;
    }
    if (predtig_serial_ == succtig_serial_)
    {
        return false;
    }
    return true;
}