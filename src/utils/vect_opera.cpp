#include <vector>
#include <numeric> // std::accumulate
#include <cmath>   // sqrt

const double CalcVectMean(const std::vector<float> &x)
{
    return (std::accumulate(x.cbegin(), x.cend(), 0.0) / x.size());
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

//#include <iostream>
const double CalcPearsonCorr_old(const std::vector<float> &x, const std::vector<float> &y)
{
    //std::cout << "vec size: " << x.size() << std::endl;
    const double mean_x = CalcVectMean(x), mean_y = CalcVectMean(y);
    double prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        prod_sum += ((x[i] - mean_x) * (y[i] - mean_y));
        t1_sqsum += ((x[i] - mean_x) * (x[i] - mean_x));
        t2_sqsum += ((y[i] - mean_y) * (y[i] - mean_y));
    }
    if (t1_sqsum == 0 || t2_sqsum == 0)
    {
        return 0;
    }
    else
    {
        return (prod_sum / sqrt(t1_sqsum * t2_sqsum));
    }
}

const double CalcPearsonCorr(const std::vector<float> &x, const std::vector<float> &y)
{
    // Compute prerequire for pearson
    double sum_x(0), sum_y(0), sum_x2(0), sum_y2(0), sum_xy(0);
    for (size_t i(0); i < x.size(); ++i)
    {
        double xi = static_cast<double>(x[i]);
        double yi = static_cast<double>(y[i]);

        sum_x += xi;
        sum_x2 += xi * xi;
        sum_y += yi;
        sum_y2 += yi * yi;
        sum_xy += xi * yi;
    }
    
    // Compute Pearson
    double prod_sum(0), t1_sqsum(0), t2_sqsum(0);
    prod_sum = sum_xy - sum_x * sum_y / x.size();
    t1_sqsum = sum_x2 - sum_x * sum_x / x.size();
    t2_sqsum = sum_y2 - sum_y * sum_y / x.size();

    // Conclude
    if (t1_sqsum == 0 || t2_sqsum == 0)
    {
        return 0;
    }
    else
    {
        return (prod_sum / sqrt(t1_sqsum * t2_sqsum));
    }
}

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y)
{
    return (0.5 * (1 - CalcPearsonCorr(x, y)));
}

const double CalcSpearmanCorr(const std::vector<float> &x, const std::vector<float> &y)
{
    static std::vector<float> x_rk, y_rk;
    CalcVectRank(x_rk, x);
    CalcVectRank(y_rk, y);
    const double spearman_corr = CalcPearsonCorr(x_rk, y_rk);
    x_rk.clear();
    y_rk.clear();
    return spearman_corr;
}

const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y)
{
    return (0.5 * (1 - CalcSpearmanCorr(x, y)));
}

const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y)
{
    const size_t nb_sample(x.size());
    double ctrst = 0.0;
    for (size_t i(0); i < nb_sample; ++i)
    {
        if (x[i] != y[i]) // if x[i] == y[i], then ctrst += 0.0
        {
            ctrst += fabs(static_cast<double>(x[i] - y[i]) / (x[i] + y[i]));
        }
    }
    return (ctrst / nb_sample);
}
