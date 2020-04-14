#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

#include <vector>
#include <cmath>

#include "mlpack/core/cv/k_fold_cv.hpp"
#include "mlpack/core/cv/metrics/f1.hpp"
#include "mlpack/methods/naive_bayes/naive_bayes_classifier.hpp"
#include "mlpack/methods/logistic_regression/logistic_regression.hpp"
#include "armadillo"
#include "boost/math/distributions/students_t.hpp"

typedef double (*evalfuncptr_t)(const arma::mat &, const arma::Row<size_t> &, size_t, size_t);

double nb_evaluation(const arma::mat &sample_counts,
                     const arma::Row<size_t> &sample_labels,
                     const size_t fold_num,
                     const size_t class_num)
{
    mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>,
                        mlpack::cv::F1<mlpack::cv::Binary>>
        eval_data(fold_num, sample_counts, sample_labels, class_num);
    double score = eval_data.Evaluate();
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double lr_evaluation(const arma::mat &sample_counts,
                     const arma::Row<size_t> &sample_labels,
                     const size_t fold_num,
                     const size_t class_num)
{
    mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>,
                        mlpack::cv::F1<mlpack::cv::Binary>>
        eval_data(fold_num, sample_counts, sample_labels);
    double score = eval_data.Evaluate();
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double mean_evaluation(const arma::mat &sample_counts,
                       const arma::Row<size_t> &sample_labels,
                       const size_t fold_num,
                       const size_t class_num)
{
    double score = arma::mean(arma::conv_to<arma::vec>::from(sample_counts));
    return score;
}

double median_evaluation(const arma::mat &sample_counts,
                         const arma::Row<size_t> &sample_labels,
                         const size_t fold_num,
                         const size_t class_num)
{
    double score = arma::median(arma::conv_to<arma::vec>::from(sample_counts));
    return score;
}

double sd_evaluation(const arma::mat &sample_counts,
                     const arma::Row<size_t> &sample_labels,
                     const size_t fold_num,
                     const size_t class_num)
{
    double score = arma::stddev(arma::conv_to<arma::vec>::from(sample_counts), 0);
    return score;
}

double rsd_evaluation(const arma::mat &sample_counts,
                      const arma::Row<size_t> &sample_labels,
                      const size_t fold_num,
                      const size_t class_num)
{
    double sd = sd_evaluation(sample_counts, sample_labels, fold_num, class_num),
           mean = mean_evaluation(sample_counts, sample_labels, fold_num, class_num);
    double score = sd / mean;
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double mc_evaluation(const arma::mat &sample_counts,
                     const arma::Row<size_t> &sample_labels,
                     const size_t fold_num,
                     const size_t class_num)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double mean1 = mean_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           mean2 = mean_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double score = (mean2 - mean1) / (mean1 + mean2);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double rsdc_evaluation(const arma::mat &sample_counts,
                       const arma::Row<size_t> &sample_labels,
                       const size_t fold_num,
                       const size_t class_num)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double rsd1 = rsd_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           rsd2 = rsd_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double score = (rsd2 - rsd1) / (rsd1 + rsd2);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double ttest_evaluation(const arma::mat &sample_counts,
                        const arma::Row<size_t> &sample_labels,
                        const size_t fold_num,
                        const size_t class_num) // mean, std, slightly different from dekupl ttestFilter, but same with R
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    for (size_t i = 0; i < cond1_counts.size(); ++i)
    {
        cond1_counts(i, 0) = log(cond1_counts(i, 0) + 1);
    }
    for (size_t i = 0; i < cond2_counts.size(); ++i)
    {
        cond2_counts(i, 0) = log(cond2_counts(i, 0) + 1);
    }
    size_t cond1_num = cond1_counts.size(), cond2_num = cond2_counts.size();
    double cond1_mean = mean_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_mean = mean_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double cond1_sd = sd_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_sd = sd_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double pvalue;
    if (cond1_sd == 0 && cond2_sd == 0)
    {
        pvalue = 1;
    }
    else
    {
        double t1, t2, df, t_stat;
        t1 = cond1_sd * cond1_sd / cond1_num;
        t2 = cond2_sd * cond2_sd / cond2_num;
        df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (cond1_num - 1) + t2 * t2 / (cond2_num - 1));
        t_stat = (cond1_mean - cond2_mean) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        pvalue = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
    }
    return pvalue;
}

double es_evaluation(const arma::mat &sample_counts,
                     const arma::Row<size_t> &sample_labels,
                     const size_t fold_num,
                     const size_t class_num)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_mean = mean_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_mean = mean_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double cond1_sd = sd_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_sd = sd_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double score = (cond2_mean - cond1_mean) / (cond1_sd + cond2_sd);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double lfcmean_evaluation(const arma::mat &sample_counts,
                          const arma::Row<size_t> &sample_labels,
                          const size_t fold_num,
                          const size_t class_num)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_mean = mean_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_mean = mean_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double score = log2(cond2_mean / cond1_mean);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

double lfcmedian_evaluation(const arma::mat &sample_counts,
                            const arma::Row<size_t> &sample_labels,
                            const size_t fold_num,
                            const size_t class_num)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0));
    arma::mat cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_median = median_evaluation(cond1_counts, sample_labels, fold_num, class_num),
           cond2_median = median_evaluation(cond2_counts, sample_labels, fold_num, class_num);
    double score = log2(cond2_median / cond1_median);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

#endif //KMEREVALUATE_EVALMETHODS_H