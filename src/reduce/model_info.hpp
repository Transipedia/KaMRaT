#ifndef KAMRAT_REDUCE_MODELINFO_HPP
#define KAMRAT_REDUCE_MODELINFO_HPP

#include <cstddef>
#include <string>
#include <iostream>

class ModelInfo
{
public:
    ModelInfo(const std::string &seq, size_t count_disk_pos, double model_score);
    const std::string &GetSeq() const;
    const size_t GetCountDiskPos() const;
    const double GetModelScore() const;
    void ScaleScore(double fact, double lower_lim, double upper_lim);

private:
    std::string seq_;
    size_t count_disk_pos_;
    double model_score_;
};

#endif //KAMRAT_REDUCE_MODELINFO_HPP