#include "feature_elem.hpp"

FeatureElem::FeatureElem(const size_t tab_serial) noexcept
    : tab_serial_(tab_serial)
{
}

const size_t FeatureElem::GetTabSerial() const noexcept
{
    return tab_serial_;
}

const double FeatureElem::GetScore() const noexcept
{
    return score_;
}

const void FeatureElem::EvalFeature(const std::vector<size_t> label_vect, const std::vector<float> &norm_count_vect)
{
    
}

const void FeatureElem::AdjustScore(const double factor, const double lower_lim, const double upper_lim)
{
    score_ *= factor;
    score_ = (score_ < lower_lim) ? lower_lim : score_;
    score_ = (score_ > upper_lim) ? upper_lim : score_;
}
