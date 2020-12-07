#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <memory>

#include "scorer.hpp"

class FeatureElem
{
public:
    FeatureElem(const size_t uniq_code, const std::unique_ptr<Scorer> &scorer);
    FeatureElem(const size_t uniq_code, const float score, const std::unique_ptr<Scorer> &scorer);
    const size_t GetUniqCode() const;
    const float GetScore() const;
    void ScaleScore(float fact, float lower_lim, float upper_lim);
    const std::vector<float> &GetCondiMeans() const;

private:
    size_t uniq_code_;
    float score_;
    std::vector<float> condi_means_;
};

using featureVect_t = std::vector<FeatureElem>;

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP