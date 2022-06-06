#include <vector>
#include <numeric> // std::accumulate
#include <cmath>   // sqrt
#include <algorithm> // sort

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


const void dual_mean_stddev(const std::vector<float> &vect, const std::vector<size_t> &categ, double &mean1, double &mean2, double &stddev1, double &stddev2, size_t &nb_catg1, size_t &nb_catg2) {
    // Compute sums for mean & stddev
    double sum1(0), sq_sum1(0), num1(0);
    double sum2(0), sq_sum2(0), num2(0);

    int idx(0);
    for (const float &x : vect) {
        double val = static_cast<double>(x);

        if (categ[idx++] == 0) {
            sum1 += val;
            sq_sum1 += val * val;
            num1 += 1;
        } else {
            sum2 += val;
            sq_sum2 += val * val;
            num2 += 1;
        }
    }

    double size1 = static_cast<double>(num1);
    double size2 = static_cast<double>(num2);
    // Compute mean
    mean1 = sum1/size1;
    mean2 = sum2/size2;
    // compute stddev
    stddev1 = sq_sum1 + (size1 * mean1 * mean1) - (2 * mean1 * sum1);
    stddev1 = sqrt(stddev1 / (size1 - 1));
    stddev2 = sq_sum2 + (size2 * mean2 * mean2) - (2 * mean2 * sum2);
    stddev2 = sqrt(stddev2 / (size2 - 1));

    nb_catg1 = num1;
    nb_catg2 = num2;
}


const void mean_stddev_min(const std::vector<float> &vect, double &mean, double &stddev, double &min) {
    // Compute sums for mean & stddev
    double sum(0), sq_sum(0);
    min = vect[0];

    for (const float &x : vect) {
        double val = static_cast<double>(x);
        sum += val;
        sq_sum += val * val;
        if (val < min) min = val;
    }

    double size = static_cast<double>(vect.size());
    // Compute mean
    mean = sum/size;
    // compute stddev
    stddev = sq_sum + (size * mean * mean) - (2 * mean * sum);
    stddev = sqrt(stddev / (size - 1));
}


const void mean_stddev(const std::vector<float> &vect, double &mean, double &stddev) {
    // Compute sums for mean & stddev
    double sum(0), sq_sum(0);

    for (const float &x : vect) {
        double val = static_cast<double>(x);
        sum += val;
        sq_sum += val * val;
    }

    double size = static_cast<double>(vect.size());
    // Compute mean
    mean = sum/size;
    // compute stddev
    stddev = sq_sum + (size * mean * mean) - (2 * mean * sum);
    stddev = sqrt(stddev / (size - 1));
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


const void getOrder(const std::vector<float> &vec, std::vector<uint> &order) {
    // Resize and init the order
    if (order.size() != vec.size())
        order.resize(vec.size());
    for (uint i=0 ; i<vec.size() ; i++)
        order[i] = i;

    // Sort order regarding vec
    std::sort(order.begin(), order.end(), [&vec](uint a, uint b){return (vec[a]<vec[b]);});
}

const void orderToRank(const std::vector<float> &vec, const std::vector<uint> & order, std::vector<float> &rank) {
    uint prev_idx = 0;
    float prev_value = vec[order[0]];
    rank.clear();
    rank.resize(vec.size());

    for (uint i=1 ; i<vec.size() ; i++) {
        if (vec[order[i]] > prev_value) {
            // Compute rank
            uint nb_similar = i - prev_idx;
            float current_rank = static_cast<float>(prev_idx) + static_cast<float>(nb_similar - 1) / 2;

            // Register ranks
            for (uint j=prev_idx ; j<i ; j++) {
                rank[order[j]] = current_rank;
            }

            // Update variables
            prev_idx = i;
            prev_value = vec[order[i]];
        }
        // Else, values are equal
    }

    // Add the last values
    // Compute rank
    uint nb_similar = vec.size() - prev_idx;
    float current_rank = static_cast<float>(prev_idx) + static_cast<float>(nb_similar - 1) / 2;

    // Register ranks
    for (uint j=prev_idx ; j<vec.size() ; j++) {
        rank[order[j]] = current_rank;
    }
}


const double CalcSpearmanCorr(const std::vector<float> &x, const std::vector<float> &y)
{
    // Get the order in both x and y vectors
    static std::vector<uint> x_order, y_order;
    getOrder(x, x_order);
    getOrder(y, y_order);
    // Transform the orders into floats
    static std::vector<float> x_rank, y_rank;
    if (x_rank.size() != x_order.size()) {
        x_rank.resize(x_order.size());
        y_rank.resize(y_order.size());
    }
    orderToRank(x, x_order, x_rank);
    orderToRank(y, y_order, y_rank);
    // Compute correlations
    const double spearman_corr = CalcPearsonCorr(x_rank, y_rank);
    return spearman_corr;
}

const double CalcSpearmanCorr_old(const std::vector<float> &x, const std::vector<float> &y)
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
