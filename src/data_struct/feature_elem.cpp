#include "feature_elem.hpp"

FeatureElem::FeatureElem(const size_t index_pos, const std::unique_ptr<Scorer> &scorer)
    : uniq_code_(index_pos)
{
    scorer->CalcCondiMeans(condi_means_);
    score_ = scorer->EvaluateScore();
}

FeatureElem::FeatureElem(const size_t uniq_code, const float score, const std::unique_ptr<Scorer> &scorer)
    : uniq_code_(uniq_code), score_(score)
{
    scorer->CalcCondiMeans(condi_means_);
}

const uint64_t FeatureElem::GetUniqCode() const
{
    return uniq_code_;
}

const float FeatureElem::GetScore() const
{
    return score_;
}

void FeatureElem::ScaleScore(float fact, float lower_lim, float upper_lim)
{
    score_ *= fact;
    score_ = (score_ < lower_lim) ? lower_lim : score_;
    score_ = (score_ > upper_lim) ? upper_lim : score_;
}

const std::vector<float> &FeatureElem::GetCondiMeans() const
{
    return condi_means_;
}