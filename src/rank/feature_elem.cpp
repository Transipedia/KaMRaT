#include <cmath>
#include <fstream>

#include "feature_elem.hpp"

FeatureElem::FeatureElem(double rep_value, std::vector<float> &count_vect, const std::string &value_str, std::ofstream &idx_file)
    : idx_pos_(static_cast<size_t>(idx_file.tellp())), score_(rep_value)
{
    idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(float)); // first index counts columns
    idx_file << value_str << std::endl;                                                          // then index non-count columns
}

const void FeatureElem::SetScore(const double score)
{
    score_ = score;
}

const double FeatureElem::GetScore() const
{
    return score_;
}

const void FeatureElem::AdjustScore(const double factor, const double lower_lim, const double upper_lim)
{
    score_ *= factor;
    score_ = (score_ < lower_lim) ? lower_lim : score_;
    score_ = (score_ > upper_lim) ? upper_lim : score_;
}

const void FeatureElem::ReserveCondiStats(const size_t nb_class)
{
    condi_mean_.reserve(nb_class);
    condi_stddev_.reserve(nb_class);
}

const void FeatureElem::AddCondiStats(const double mean, const double stddev)
{
    condi_mean_.push_back(mean);
    condi_stddev_.push_back(stddev);
}

const double FeatureElem::GetCondiStats(const size_t i_condi, const std::string &&stats_name) const
{
    if (stats_name == "mean")
    {
        return condi_mean_[i_condi];
    }
    else if (stats_name == "stddev")
    {
        return condi_stddev_[i_condi];
    }
    else
    {
        return std::nan("");
    }
}

const void FeatureElem::RetrieveCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count) const
{
    count_vect.resize(nb_count);
    idx_file.seekg(idx_pos_);
    idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(float));
}
