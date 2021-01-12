#ifndef KAMRAT_RANK_FEATUREELEM_HPP
#define KAMRAT_RANK_FEATUREELEM_HPP

#include <vector>
#include <tuple>
#include <fstream>

class FeatureElem
{
public:
    FeatureElem(double rep_value, std::vector<float> &count_vect, const std::string &value_str, std::ofstream &idx_file);

    const void SetScore(double score);
    const double GetScore() const;
    const void AdjustScore(double factor, double lower_lim, double upper_lim);

    const void ReserveCondiStats(size_t nb_class);
    const void AddCondiStats(double mean, double stddev);
    const std::vector<double> &GetCondiMeanVect() const;
    const std::vector<double> &GetCondiStddevVect() const;
    const double GetCondiMeanAt(size_t i_smp) const;
    const double GetCondiStddevAt(size_t i_smp) const;
    const void RetrieveCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, size_t nb_count) const;
    const void RetrieveCountVectValueStr(std::vector<float> &count_vect, std::string &value_str, std::ifstream &idx_file, const size_t nb_count) const;

private:
    size_t idx_pos_;                                // serial in count table
    double score_;                                  // feature's score
    std::vector<double> condi_mean_, condi_stddev_; // condition's mean, median, stddev
};

using featureVect_t = std::vector<FeatureElem>;

#endif //KAMRAT_RANK_FEATUREELEM_HPP