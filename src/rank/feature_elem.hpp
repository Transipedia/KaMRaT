#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <memory>
#include <fstream>

#include "tab_elem.hpp"

class FeatureElem
{
public:
    FeatureElem(size_t tab_serial) noexcept;

    const size_t GetTabSerial() const noexcept;
    const double GetScore() const noexcept;
    const void EvalFeature(const std::vector<size_t> label_vect, const TabElem &tab_elem, const std::vector<float> &nf_vect,
                           const std::string &score_method, std::ifstream &idx_file);
    const void AdjustScore(double factor, double lower_lim, double upper_lim);

private:
    size_t tab_serial_;                                  // serial in count table
    double score_;                                       // feature's score
    std::vector<std::pair<float, float>> condi_mean_sd_; // mean and standard deviation of each condition
};

using featureVect_t = std::vector<FeatureElem>;

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP