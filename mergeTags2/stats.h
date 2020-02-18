//
// Created by haoliang.xue on 6/28/19.
//

#ifndef STATS_H
#define STATS_H

#include <math.h>
#include <assert.h>

#include "alglib/statistics.h"
#include "alglib/interpolation.h"

#define MET_NONE "none"
#define MET_PEARSON "pearson"
#define MET_SPEARMAN "spearman"
#define MET_CONTRAST "contrast"
#define MIN_PEARSON 0.61
#define MIN_SPEARMAN 0.56
#define MAX_CONTRAST 0.25

static inline float absolute(float x)
{
    return x < 0 ? -x : x;
}

static inline int sign(float x)
{
    if (x < 0)
        return -1;
    else if (x > 0)
        return 1;
    else
        return 0;
}

static inline float mean(std::vector<float> x)
{
    float sum = 0;
    for (float xi : x)
    {
        sum += xi;
    }
    return (sum / x.size());
}

static inline float cov(std::vector<float> x, std::vector<float> y)
{
    float mean_x = mean(x), mean_y = mean(y), cov = 0;
    size_t n = x.size();
    for (size_t i = 0; i < n; i++)
    {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return (cov / n);
}

static inline float var(std::vector<float> x)
{
    float mean_x = mean(x), var_x = 0;
    for (float xi : x)
    {
        var_x += (xi - mean_x) * (xi - mean_x);
    }
    return var_x;
}

static inline float sd(std::vector<float> x)
{
    return sqrt(var(x) / x.size());
}

static inline float calc_pearson_corr(const float *x, const float *y, const size_t n_sample)
{
    double *xx = (double *)calloc(n_sample, sizeof(double)), *yy = (double *)calloc(n_sample, sizeof(double));
    for (size_t i = 0; i < n_sample; i ++) {
        xx[i] = x[i];
        yy[i] = y[i];
    }
    alglib::real_1d_array x1, y1;
    x1.setcontent(n_sample, xx);
    y1.setcontent(n_sample, yy);
    free(xx);
    free(yy);
    return alglib::pearsoncorr2(x1, y1);
}

static inline float calc_spearman_corr(const float *x, const float *y, const size_t n_sample)
{
    double *xx = (double *)calloc(n_sample, sizeof(double)), *yy = (double *)calloc(n_sample, sizeof(double));
    for (size_t i = 0; i < n_sample; i ++) {
        xx[i] = x[i];
        yy[i] = y[i];
    }
    alglib::real_1d_array x1, y1;
    x1.setcontent(n_sample, xx);
    y1.setcontent(n_sample, yy);
    free(xx);
    free(yy);
    return alglib::spearmancorr2(x1, y1);
}

static inline float calc_proportionality(const float *x, const float *y, const size_t n_sample)
{
    std::vector<float> log_x(n_sample), log_y(n_sample), log_d(n_sample);
    for (size_t i = 0; i < n_sample; i++)
    {
        log_x[i] = log(x[i] + 1E-3);
        log_y[i] = log(y[i] + 1E-3);
        log_d[i] = log_x[i] - log_y[i];
    }
    assert(var(log_x) != 0);
    return (var(log_d) / var(log_x));
}

static inline float calc_mean_contrast(const float *x, const float *y, const size_t n_sample)
{
    std::vector<float> ctrst(n_sample);
    for (size_t i = 0; i < n_sample; i++)
    {
        if (x[i] + y[i] > 0) {
            ctrst[i] = absolute((x[i] - y[i]) / (x[i] + y[i]));
        } else {
            ctrst[i] = 0;
        }
    }
    return (mean(ctrst));
}

#endif //STATS_H
