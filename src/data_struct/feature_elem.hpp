#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <memory>

#include "scorer.hpp"
#include "tab_elem.hpp"

class FeatureElem : public TabElem
{
public:
    FeatureElem(std::istringstream &line_conv, std::vector<float> &count_vect, std::ofstream &idx_file, const TabHeader &tab_header);

    const void EvalFeatureElem(const std::unique_ptr<Scorer> &scorer);
    void ScaleValue(float fact, float lower_lim, float upper_lim);
    const std::vector<float> &GetNormCondiMeans() const;

private:
    std::vector<float> norm_condi_means_;
};

using featureVect_t = std::vector<FeatureElem>;

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP