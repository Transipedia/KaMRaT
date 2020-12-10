#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <memory>

#include "scorer.hpp"
#include "tab_elem.hpp"

class FeatureElem
{
public:
    const void MakeFeatureElem(const std::unique_ptr<Scorer> &scorer, float score_from_tab);
    const float GetScore() const;
    void ScaleScore(float fact, float lower_lim, float upper_lim);
    const std::vector<float> &GetCondiMeans() const;

private:
    float score_;
    std::vector<float> condi_means_;
};

using featureVect_t = std::vector<std::pair<size_t, FeatureElem>>;

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP