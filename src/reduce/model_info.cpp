#include "model_info.hpp"

ModelInfo::ModelInfo(const std::string &seq, size_t count_disk_pos, double model_score)
    : seq_(seq), count_disk_pos_(count_disk_pos), model_score_(model_score)
{
}

const std::string &ModelInfo::GetSeq() const
{
    return seq_;
}
const size_t ModelInfo::GetCountDiskPos() const
{
    return count_disk_pos_;
}

const double ModelInfo::GetModelScore() const
{
    return model_score_;
}

void ModelInfo::ScaleScore(double fact, double lower_lim, double upper_lim)
{
    model_score_ *= fact;
    model_score_ = (model_score_ < lower_lim) ? lower_lim : model_score_;
    model_score_ = (model_score_ > upper_lim) ? upper_lim : model_score_;
}