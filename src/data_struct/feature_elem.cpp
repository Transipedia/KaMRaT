#include "feature_elem.hpp"

FeatureElem::FeatureElem(std::vector<float> &count_vect, std::string &non_count_str, const float value, std::ofstream &idx_file)
    : TabElem(count_vect, non_count_str, value, idx_file)
{
}

const void FeatureElem::EvalFeatureElem(const std::unique_ptr<Scorer> &scorer)
{
    scorer->CalcNormCondiMeans(norm_condi_means_);
    if (scorer->GetScoreMethod() != "user")
    {
        scorer->TransformCounts();
        value_ = scorer->EvaluateScore();
    }
}

void FeatureElem::ScaleValue(float fact, float lower_lim, float upper_lim)
{
    value_ *= fact;
    value_ = (value_ < lower_lim) ? lower_lim : value_;
    value_ = (value_ > upper_lim) ? upper_lim : value_;
}

const std::vector<float> &FeatureElem::GetNormCondiMeans() const
{
    return norm_condi_means_;
}