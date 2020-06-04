#ifndef KAMRAT_UTILS_STATISTICS_HPP
#define KAMRAT_UTILS_STATISTICS_HPP

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <map>
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
inline double CalcVectMean(const std::vector<countT> &x)
{
    double sum = std::accumulate(x.cbegin(), x.cend(), 0.0);
    return (sum / x.size());
}

template <typename countT>
inline void CalcVectRank(std::vector<float> &x_rk, const std::vector<countT> &x)
{
    x_rk.resize(x.size());
    std::multimap<float, size_t> x_pos; // utilize map's automatic orderign
    for (size_t i(0); i < x.size(); ++i)
    {
        x_pos.insert({x.at(i), i});
    }
    int rk = 1;
    for (auto iter_x_pos = x_pos.cbegin(); iter_x_pos != x_pos.cend(); iter_x_pos = x_pos.upper_bound(iter_x_pos->first))
    {
        size_t c = x_pos.count(iter_x_pos->first);
        for (auto iter_range = x_pos.equal_range(iter_x_pos->first).first; iter_range != x_pos.equal_range(iter_x_pos->first).second; ++iter_range)
        {
            x_rk.at(iter_range->second) = (2 * rk + c - 1) / 2.0;
        }
        rk += c;
    }
}

template <typename countT>
inline double CalcPearsonCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    double mean_x = CalcVectMean(x), mean_y = CalcVectMean(y), prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        double t1 = x.at(i) - mean_x, t2 = y.at(i) - mean_y;
        prod_sum += (t1 * t2);
        t1_sqsum += (t1 * t1);
        t2_sqsum += (t2 * t2);
    }
    return (prod_sum / sqrt(t1_sqsum * t2_sqsum));
}

template <typename countT>
inline double CalcSpearmanCorrelation(const std::vector<countT> &x, const std::vector<countT> &y)
{
    std::vector<float> x_rk, y_rk;
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
    else
    {
        return -1;
    }
}

#endif //KAMRAT_UTILS_STATISTICS_HPP
