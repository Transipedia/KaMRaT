#ifndef KAMRAT_RANK_FEATUREELEM_HPP
#define KAMRAT_RANK_FEATUREELEM_HPP

#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include <iostream>

class FeatureElem
{
public:
    FeatureElem(const std::string &feature, size_t pos);
    FeatureElem(const std::string &feature, const std::vector<size_t> &mem_pos_vect);

    const std::string &GetFeature() const;
    const size_t GetRepPos() const;
    const size_t GetNbMemPos() const;
    const double GetScore() const;
    void SetScore(double score);
    const double AdjustScore(double factor, double lower_lim, double upper_lim);
    const std::vector<float> &EstimateCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, size_t nb_smp, const std::string &count_mode) const;

    static double AdjustScore(double score, const double factor, const double lower_lim, const double upper_lim)
    {
        auto save = score;
        score *= factor;
        score = (score < lower_lim) ? lower_lim : score;
        score = (score > upper_lim) ? upper_lim : score;

        if (score == 0)
            std::cout << save << " -> 0" << std::endl;
        return score;
    }

private:
    const std::string feature_;              // feature name or contig sequence
    const std::vector<size_t> mem_pos_vect_; // member k-mer position vector, the first is rep-pos (for k-mers and general features)
    double score_;                           // feature score
};

#endif //KAMRAT_RANK_FEATUREELEM_HPP