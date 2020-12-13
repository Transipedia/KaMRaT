#ifndef KAMRAT_UTILS_VECOPERATION_HPP
#define KAMRAT_UTILS_VECOPERATION_HPP

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <set>
#include <iterator>

template <typename countT>
inline float CalcVectMean(const std::vector<countT> &x)
{
    float sum = std::accumulate(x.cbegin(), x.cend(), 0.0);
    return (sum / x.size());
}

template <typename countT>
inline float CalcVectMedian(const std::vector<countT> &x)
{
    static std::vector<countT> x_tmp;
    x_tmp = x;
    std::sort(x_tmp.begin(), x_tmp.end());
    const size_t nb_count = x.size();
    float median = ((nb_count % 2) ? x_tmp[nb_count / 2] : (0.5 * (x_tmp[nb_count / 2] + x_tmp[nb_count / 2 + 1])));
    x_tmp.clear();
    return median;
}

template <typename countT>
inline void CalcVectRank(std::vector<float> &x_rk, const std::vector<countT> &x)
{
    static std::vector<size_t> r, s; // r for rank number, s for same number
    size_t n = x.size();
    r.resize(n, 1);
    s.resize(n, 1);
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
    r.clear();
    s.clear();
}

template <typename countT>
inline float CalcPearsonDist(const std::vector<countT> &x, const std::vector<countT> &y)
{
    float mean_x = CalcVectMean(x), mean_y = CalcVectMean(y), prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        prod_sum += ((x[i] - mean_x) * (y[i] - mean_y));
        t1_sqsum += ((x[i] - mean_x) * (x[i] - mean_x));
        t2_sqsum += ((y[i] - mean_y) * (y[i] - mean_y));
    }
    return (0.5 * (1 - (prod_sum / sqrt(t1_sqsum * t2_sqsum))));
}

template <typename countT>
inline float CalcSpearmanDist(const std::vector<countT> &x, const std::vector<countT> &y)
{
    static std::vector<float> x_rk, y_rk;
    CalcVectRank(x_rk, x);
    CalcVectRank(y_rk, y);
    x_rk.clear();
    y_rk.clear();
    return CalcPearsonDist(x_rk, y_rk);
}

template <typename countT>
inline float CalcMACDist(const std::vector<countT> &x, const std::vector<countT> &y)
{
    static std::vector<float> ctrst;
    const size_t nb_sample(x.size());
    ctrst.resize(nb_sample);
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
    ctrst.clear();
    return CalcVectMean(ctrst);
}

template <typename countT>
inline float CalcXDist(const std::vector<countT> &x, const std::vector<countT> &y, const std::string &eval_method)
{
    if (x.size() != y.size())
    {
        throw std::domain_error("two vectors have different size");
    }
    if (eval_method == "mac")
    {
        return CalcMACDist(x, y);
    }
    else if (eval_method == "pearson")
    {
        return CalcPearsonDist(x, y);
    }
    else if (eval_method == "spearman")
    {
        return CalcSpearmanDist(x, y);
    }
    else if (eval_method == "none")
    {
        return 0;
    }
    else
    {
        throw std::invalid_argument("unknown intervention/evaluation method: " + eval_method);
    }
}

template <typename countT>
const std::vector<countT> &CalcVectsRes(std::vector<countT> &vect_res, const std::vector<std::vector<countT>> &vect_vect, const std::string &method)
{
    if (method == "mean")
    {
        for (const auto &v : vect_vect)
        {
            vect_res.emplace_back(CalcVectMean(v));
        }
    }
    else if (method == "median")
    {
        for (const auto &v : vect_vect)
        {
            vect_res.emplace_back(CalcVectMedian(v));
        }
    }
    else // should not happen, for debug
    {
        throw std::domain_error("unknown method: " + method);
    }
    return vect_res;
}

#endif //KAMRAT_UTILS_VECOPERATION_HPP