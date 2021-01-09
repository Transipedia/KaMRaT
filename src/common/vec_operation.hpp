#ifndef KAMRAT_COMMON_VECOPERATION_HPP
#define KAMRAT_COMMON_VECOPERATION_HPP

#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

/* ======================================= *\
 * Used by kamratMerge and contigEvaluator *
\* ======================================= */

const double CalcVectMean(const std::vector<float> &x)
{
    return (std::accumulate(x.cbegin(), x.cend(), 0.0) / x.size());
}

const double CalcVectMedian(const std::vector<float> &x)
{
    static std::vector<float> x_tmp;
    x_tmp = x;
    size_t n_elem = x_tmp.size();
    auto mid_elem = x_tmp.begin() + n_elem / 2;
    std::nth_element(x_tmp.begin(), mid_elem, x_tmp.end());
    if (n_elem % 2 != 0)
    {
        return *mid_elem;
    }
    else
    {
        return 0.5 * (*mid_elem + *(std::max_element(x_tmp.begin(), mid_elem)));
    }
}

void CalcVectRank(std::vector<float> &x_rk, const std::vector<float> &x)
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

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y)
{
    const double mean_x = CalcVectMean(x), mean_y = CalcVectMean(y);
    double prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        prod_sum += ((x[i] - mean_x) * (y[i] - mean_y));
        t1_sqsum += ((x[i] - mean_x) * (x[i] - mean_x));
        t2_sqsum += ((y[i] - mean_y) * (y[i] - mean_y));
    }
    return (0.5 * (1 - (prod_sum / sqrt(t1_sqsum * t2_sqsum))));
}

const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y)
{
    static std::vector<float> x_rk, y_rk;
    CalcVectRank(x_rk, x);
    CalcVectRank(y_rk, y);
    x_rk.clear();
    y_rk.clear();
    return CalcPearsonDist(x_rk, y_rk);
}

const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y)
{
    const size_t nb_sample(x.size());
    double ctrst = 0.0;
    for (size_t i(0); i < nb_sample; ++i)
    {
        if (x[i] != y[i]) // if x[i] == y[i], then ctrst += 0.0
        {
            ctrst += fabs(static_cast<double>(x[i] - y[i])) / (x[i] + y[i]);
        }
    }
    return (ctrst / nb_sample);
}

#endif //KAMRAT_COMMON_VECOPERATION_HPP