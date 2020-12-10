#include "feature_elem.hpp"

const void FeatureElem::MakeFeatureElem(const std::unique_ptr<Scorer> &scorer, const float score_from_tab)
{
    if (scorer->GetScoreMethod() != "user")
    {
        scorer->TransformCounts();
        score_ = scorer->EvaluateScore();
    }
    else
    {
        score_ = score_from_tab;
    }
}

const uint64_t FeatureElem::GetIdxPos() const
{
    return idx_pos_;
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