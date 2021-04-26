#ifndef KAMRAT_RANK_FEATUREELEM_HPP
#define KAMRAT_RANK_FEATUREELEM_HPP

#include <string>
#include <vector>
#include <memory>
#include <fstream>

class FeatureElem
{
public:
    FeatureElem(const std::string &&feature, size_t pos);
    FeatureElem(const std::string &&feature, const std::vector<size_t> &&mem_pos_vect);

    const double GetScore() const;
    const double AdjustScore(double factor, double lower_lim, double upper_lim);
    const std::vector<float> &EstimateCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, size_t nb_smp, const std::string &&count_mode) const;

private:
    const std::string feature_;              // feature name or contig sequence
    const std::vector<size_t> mem_pos_vect_; // member k-mer position vector, the first is rep-pos (for k-mers and general features)
    double score_;                           // feature score
};

using featureVect_t = std::vector<std::unique_ptr<FeatureElem>>;

#endif //KAMRAT_RANK_FEATUREELEM_HPP