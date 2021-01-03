#include <cmath>

#include "feature_elem.hpp"

FeatureElem::FeatureElem(const size_t idx_pos, const double rep_value)
    : idx_pos_(idx_pos), score_(rep_value)
{
}

const void FeatureElem::OpenIdxFile(const std::string &&idx_path)
{
    idx_file_.open(idx_path);
    if (!idx_file_.is_open())
    {
        throw std::domain_error("error open index file: " + idx_path);
    }
}

const void FeatureElem::CloseIdxFile()
{
    idx_file_.close();
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

const void FeatureElem::RetrieveCountVect(std::vector<double> &count_vect, const size_t nb_count) const
{
    count_vect.resize(nb_count);
    idx_file_.seekg(idx_pos_);
    idx_file_.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(double));
}
