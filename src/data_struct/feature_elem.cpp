#include "feature_elem.hpp"

const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp
const std::vector<float> &GetMeanCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                           const std::vector<size_t> mem_pos_vect); // in utils/index_loading.cpp
const std::vector<float> &GetMedianCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                             const std::vector<size_t> mem_pos_vect); // in utils/index_loading.cpp

FeatureElem::FeatureElem(const std::string &&feature, const size_t pos)
    : feature_(feature), mem_pos_vect_(1, pos)
{
}

FeatureElem::FeatureElem(const std::string &&feature, const std::vector<size_t> &&mem_pos_vect)
    : feature_(feature), mem_pos_vect_(mem_pos_vect)
{
}

const double FeatureElem::GetScore() const
{
    return score_;
}

const double FeatureElem::AdjustScore(const double factor, const double lower_lim, const double upper_lim)
{
    score_ *= factor;
    score_ = (score_ < lower_lim) ? lower_lim : score_;
    score_ = (score_ > upper_lim) ? upper_lim : score_;
}

const std::vector<float> &FeatureElem::EstimateCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat,
                                                         const size_t nb_smp, const std::string &&count_mode) const
{
    if (mem_pos_vect_.size() == 1 || count_mode == "rep")
    {
        GetCountVect(count_vect, idx_mat, mem_pos_vect_[0], nb_smp);
    }
    else if (count_mode == "mean")
    {
        GetMeanCountVect(count_vect, idx_mat, nb_smp, mem_pos_vect_);
    }
    else // count_mode == "median"
    {
        GetMedianCountVect(count_vect, idx_mat, nb_smp, mem_pos_vect_);
    }
}
