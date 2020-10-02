#ifndef KAMRAT_UTILS_STATISTICS_HPP
#define KAMRAT_UTILS_STATISTICS_HPP

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <set>
#include <iterator>

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
inline float CalcVectMean(const std::vector<countT> &x)
{
    float sum = std::accumulate(x.cbegin(), x.cend(), 0.0);
    return (sum / x.size());
}

template <typename countT>
inline void CalcVectRank(std::vector<float> &x_rk, const std::vector<countT> &x)
{
    size_t n = x.size();
    std::vector<size_t> r(n, 1), s(n, 1); // r for rank number, s for same number
    for (size_t i(0); i < n; ++i)
    {
        for (size_t j(i + 1); j < n; ++j)
        {
            if (x[i] == x[j])
            {
                ++s[i];
                ++s[j];
            }
            else if (x[i] < x[j])
            {
                ++r[j];
            }
            else // x[i] > x[j]
            {
                ++r[i];
            }
        }
        x_rk.push_back(r[i] + 0.5 * (s[i] - 1));
    }
}

template <typename countT>
inline float CalcPearsonCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    float mean_x = CalcVectMean(x), mean_y = CalcVectMean(y), prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        float t1 = x[i] - mean_x, t2 = y[i] - mean_y;
        prod_sum += (t1 * t2);
        t1_sqsum += (t1 * t1);
        t2_sqsum += (t2 * t2);
    }
    return (prod_sum / sqrt(t1_sqsum * t2_sqsum));
}

template <typename countT>
inline float CalcSpearmanCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    std::vector<float> x_rk, y_rk;
    CalcVectRank(x_rk, x);
    CalcVectRank(y_rk, y);
    return CalcPearsonCorrelation(x_rk, y_rk);
}

template <typename countT>
inline float CalcMeanAbsoluteContrast(const std::vector<countT> &x, const std::vector<countT> &y)
{
    const size_t nb_sample(x.size());
    std::vector<float> ctrst(nb_sample);
    for (size_t i(0); i < nb_sample; ++i)
    {
        if (x[i] == y[i]) // including the case when x[i] = y[i] = 0
        {
            ctrst[i] = 0.0;
        }
        else
        {
            ctrst[i] = fabs(static_cast<float>(x[i] - y[i])) / (x[i] + y[i]);
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
    if (x.size() != y.size())
    {
        throw std::domain_error("two vectors have different size");
    }
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
    else
    {
        throw std::invalid_argument("unknown intervention/evaluation method: " + eval_method);
    }
}

#endif //KAMRAT_UTILS_STATISTICS_HPP
