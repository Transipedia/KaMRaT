#ifndef KAMRAT_UTILS_STATISTICS_HPP
#define KAMRAT_UTILS_STATISTICS_HPP

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

template <typename countT>
inline countT GetMaxInPair(countT x, countT y)
{
    return (x > y) ? x : y;
}

template <typename countT>
inline countT GetMinInPair(countT x, countT y)
{
    return (x < y) ? x : y;
}

template <typename countT>
inline double CalcVectMean(const std::vector<countT> &x)
{
    double sum = std::accumulate(x.cbegin(), x.cend(), 0.0);
    return (sum / x.size());
}

template <typename countT>
inline void CalcVectRank(std::vector<size_t> &x_rk, const std::vector<countT> &x)
{
    std::vector<size_t> x_idx(x.size());
    std::iota(x_idx.begin(), x_idx.end(), 0);
    std::stable_sort(x_idx.begin(),
                     x_idx.end(),
                     [&x](const size_t i, const size_t j) { return x.at(i) < x.at(j); });
    x_rk.resize(x.size());
    std::iota(x_rk.begin(), x_rk.end(), 0);
    std::stable_sort(x_rk.begin(),
                     x_rk.end(),
                     [&x_idx](const size_t i, const size_t j) { return x_idx.at(i) < x_idx.at(j); });
}

template <typename countT>
inline double CalcPearsonCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    double mean_x = CalcVectMean(x), mean_y = CalcVectMean(y), prod_sum(0), t1_sqsum(0), t2_sqsum(0);

    // for (countT xx : x)
    // {
    //     std::cout << xx << "\t";
    // }
    // std::cout << "=>\t" << mean_x << std::endl;

    // for (countT yy : y)
    // {
    //     std::cout << yy << "\t";
    // }
    // std::cout << "=>\t" << mean_y << std::endl;

    for (size_t i(0); i < x.size(); ++i)
    {
        double t1 = x.at(i) - mean_x, t2 = y.at(i) - mean_y;
        // std::cout << "\t" << t1 << "," << t2;
        prod_sum += (t1 * t2);
        t1_sqsum += (t1 * t1);
        t2_sqsum += (t2 * t2);
    }
    // std::cout << "\t=>\t" << prod_sum << "," << t1_sqsum << "," << t2_sqsum << "\t=>\t" << prod_sum / sqrt(t1_sqsum * t2_sqsum) << std::endl;
    return (prod_sum / sqrt(t1_sqsum * t2_sqsum));
}

template <typename countT>
inline double CalcSpearmanCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    std::vector<size_t> x_rk, y_rk;
    CalcVectRank(x_rk, x);
    CalcVectRank(y_rk, y);
    return CalcPearsonCorrelation(x_rk, y_rk);
}

template <typename countT>
inline double CalcMeanAbsoluteContrast(const std::vector<countT> &x, const std::vector<countT> &y)
{
    const size_t nb_sample(x.size());
    std::vector<double> ctrst(nb_sample);
    for (size_t i(0); i < nb_sample; ++i)
    {
        if (x.at(i) == y.at(i)) // including the case when x[i] = y[i] = 0
        {
            ctrst.at(i) = 0.0;
        }
        else
        {
            ctrst.at(i) = fabs(static_cast<double>(x.at(i) - y.at(i))) / (x.at(i) + y.at(i));
        }
    }
    return (CalcVectMean(ctrst));
}

template <typename countT>
inline float CalcPearsonDistance(const std::vector<countT> &x, const std::vector<countT> &y)
{
    return 0.5 * (1 - CalcPearsonCorrelation(x, y));
}

template <typename countT>
inline float CalcSpearmanDistance(const std::vector<countT> &x, const std::vector<countT> &y)
{
    return 0.5 * (1 - CalcSpearmanCorrelation(x, y));
}

template <typename countT>
inline float CalcMACDistance(const std::vector<countT> &x, const std::vector<countT> &y)
{
    return CalcMeanAbsoluteContrast(x, y);
}

template <typename countT>
inline float CalcDistance(const std::vector<countT> &x, const std::vector<countT> &y, const std::string &eval_method)
{
    if (eval_method == "mac")
    {
        return CalcMACDistance(x, y);
    }
    else if (eval_method == "pearson")
    {
        return CalcPearsonDistance(x, y);
    }
    else if (eval_method == "spearman")
    {
        return CalcSpearmanDistance(x, y);
    }
}

#endif //KAMRAT_UTILS_STATISTICS_HPP
