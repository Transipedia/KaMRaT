#include "feature_elem.hpp"

FeatureElem::FeatureElem(std::istringstream &line_conv, std::vector<float> &count_vect, std::ofstream &idx_file, const TabHeader &tab_header)
    : TabElem(line_conv, count_vect, idx_file, tab_header)
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

void FeatureElem::ScaleValue(double fact, double lower_lim, double upper_lim)
{
    value_ *= fact;
    value_ = (value_ < lower_lim) ? lower_lim : value_;
    value_ = (value_ > upper_lim) ? upper_lim : value_;
}

const std::vector<float> &FeatureElem::GetNormCondiMeans() const
{
    return norm_condi_means_;
}