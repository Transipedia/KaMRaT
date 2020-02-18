#include <sstream>

#include <mlpack/core/cv/k_fold_cv.hpp>

#include "kmer_model.h"

KmerModel::KmerModel(const size_t sample_num,
                     const std::string &count_line,
                     const std::string &header_line,
                     const arma::Row<size_t> &sample_labels,
                     const std::unordered_map<std::string, size_t> &sample_info,
                     size_t fold_num,
                     const evalfuncptr_t eval_func,
                     const std::string &user_method_name)
{
    std::istringstream header_conv(header_line), count_conv(count_line);
    std::string header_x;
    float count_x;

    header_conv >> header_x;
    count_conv >> kmer_seq_;
    if (eval_func)
    {
        arma::mat sample_counts(1, sample_num, arma::fill::zeros);
        size_t i = 0;
        
        for (count_conv >> count_x, header_conv >> header_x; i < sample_num && !count_conv.fail(); count_conv >> count_x, header_conv >> header_x)
        {
            if (sample_info.empty() || sample_info.find(header_x) != sample_info.cend())
            {
                sample_counts(0, i) = count_x;
                i++;
            }
        }
        score_ = eval_func(sample_counts, sample_labels, fold_num, CLASS_NUM);
    }
    else // if using user-defined method, no need to load sample counts
    {
        for (count_conv >> count_x, header_conv >> header_x; !count_conv.fail(); count_conv >> count_x, header_conv >> header_x)
        {
            if (header_x == user_method_name)
            {
                score_ = count_x;
                break;
            }
        }
        if (count_conv.fail())
        {
            std::cerr << "ERROR: no column corresponded to evaluation method " << user_method_name << "." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

void KmerModel::Print() const
{
    std::cout << kmer_seq_ << "\t"
              << score_ << std::endl;
}

const double KmerModel::GetScore() const
{
    return score_;
}

void KmerModel::ScaleScore(const double factor, const double lower_limit, const double upper_limit)
{
    score_ *= factor;
    if (score_ < lower_limit)
    {
        score_ = lower_limit;
    }
    else if (score_ > upper_limit)
    {
        score_ = upper_limit;
    }
}