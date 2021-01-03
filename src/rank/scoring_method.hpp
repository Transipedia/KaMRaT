
const double CalcVectMean(const std::vector<float> &x)
{
    return (std::accumulate(x.cbegin(), x.cend(), 0.0) / x.size());
}

const void CalcVectMeanStddev(double &x_mean, double &x_stddev, const std::vector<float> &x)
{
    size_t x_len = x.size();
    x_mean = std::accumulate(x.cbegin(), x.cend(), 0.0) / x_len;
    double sse = 0;
    for (size_t i(0); i < x_len; ++i)
    {
        sse += pow(x[i] - x_mean, 2);
    }
    x_stddev = sqrt(sse / (x_len - 1));
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