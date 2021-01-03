#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <tuple>
#include <fstream>

class FeatureElem
{
public:
    // FeatureElem(double rep_value, std::vector<double> &count_vect, const std::string &value_str, std::ofstream &idx_file);
    FeatureElem(size_t idx_pos_, double rep_value);

    const void OpenIdxFile(const std::string &&idx_path);
    const void CloseIdxFile();

    const void SetScore(double score);
    const double GetScore() const;
    const void AdjustScore(double factor, double lower_lim, double upper_lim);

    const void ReserveCondiStats(size_t nb_class);
    const void AddCondiStats(double mean, double stddev);
    const double GetCondiStats(size_t i_condi, const std::string &&stats_name) const;
    const void RetrieveCountVect(std::vector<double> &count_vect, size_t nb_count) const;

private:
    static std::ifstream idx_file_;                 // file stream for retrieving count vector from index, all instances share one (static)
    size_t idx_pos_;                                // serial in count table
    double score_;                                  // feature's score
    std::vector<double> condi_mean_, condi_stddev_; // condition's mean, median, stddev
};

using featureVect_t = std::vector<FeatureElem>;

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP