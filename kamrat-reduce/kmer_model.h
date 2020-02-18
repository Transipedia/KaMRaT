#ifndef KMEREVALUATE_KMERMODEL_H
#define KMEREVALUATE_KMERMODEL_H

#include <string>

#include <mlpack/core.hpp>

#include "utils.h"

class KmerModel
{
public:
    KmerModel(size_t sample_num,
              const std::string &count_line,
              const std::string &header_line,
              const arma::Row<size_t> &sample_labels,
              const std::unordered_map<std::string, size_t> &sample_info,
              size_t fold_num,
              const evalfuncptr_t eval_func,
              const std::string &user_method_name);
    void Print() const;
    const double GetScore() const;
    void ScaleScore(double factor, double lower_limit, double upper_limit);

private:
    std::string kmer_seq_;
    double score_ = -1;
};

#endif //KMEREVALUATE_KMERMODEL_H