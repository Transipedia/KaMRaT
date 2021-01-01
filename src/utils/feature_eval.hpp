#ifndef KAMRAT_UTILS_FEATUREEVAL_HPP
#define KAMRAT_UTILS_FEATUREEVAL_HPP

#include <vector>
#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include <mlpack/core/cv/k_fold_cv.hpp>
#include <mlpack/core/cv/metrics/f1.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/methods/linear_svm/linear_svm.hpp>
#include <mlpack/methods/linear_svm/linear_svm_function.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <armadillo>

/* ============================ Scoring Method ============================ *\
 * rsd          relative standard deviation                                  *
 * ttest        p-value of t-test, adjusted by Benjamini-Hochberg procedure  *
 * snr          signal-to-noise ratio                                        *
 * lfc:mean     log2-fold-change between condition means                     *
 * lfc:median   Log2-fold-change between condition medians                   *
 * nbc          naive Bayes classifier                                       *
 * lrc          logistic regression classifier                               *
 * svm          support vector machine                                       *
 * user:name    use representative value in table as score                   *
\* ======================================================================== */

class MLMetrics
{
public:
    template <typename MLAlgorithm, typename DataType>
    static double Evaluate(MLAlgorithm &model, const DataType &data, const arma::Row<size_t> &labels)
    {
        if (labels.max() == 1) // binary conditions
        {
            double score = mlpack::cv::F1<mlpack::cv::Binary>().Evaluate(model, data, labels);
            return (isnan(score) ? 0 : score);
        }
        else // multiple conditions
        {
            double score = mlpack::cv::F1<mlpack::cv::Micro>().Evaluate(model, data, labels);
            return score;
        }
    }
};

const double CalcRelatSDScore(const double all_mean, const double all_stddev)
{
    if (all_mean <= 1)
    {
        return all_stddev;
    }
    else
    {
        return (all_stddev / all_mean);
    }
}

const double CalcTtestScore(const size_t nb1, const size_t nb2,
                            const double mean1, const double mean2,
                            const double stddev1, const double stddev2)
{
    double pvalue;
    if (stddev1 == 0 && stddev2 == 0)
    {
        pvalue = 1;
    }
    else
    {
        double t1 = stddev1 * stddev1 / nb1,
               t2 = stddev2 * stddev2 / nb2,
               df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (nb1 - 1) + t2 * t2 / (nb2 - 1)),
               t_stat = (mean1 - mean2) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        pvalue = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
    }
    return pvalue;
}

const double CalcSNRScore(const double mean1, const double mean2,
                          const double stddev1, const double stddev2)
{
    return ((mean1 - mean2) / (stddev1 + stddev2));
}

const double CalcNBCScore(const size_t nb_fold, const size_t nb_class,
                          const std::vector<size_t> &label_vect, const std::vector<float> &norm_count_vect)
{
    const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(norm_count_vect, label_vect, nb_class);
        MLMetrics my_f1;
        return my_f1.Evaluate(nbc, norm_count_vect, label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
            score_data(nb_fold_final, norm_count_vect, label_vect, nb_class);
        return score_data.Evaluate();
    }
}

const double CalcLRCScore(const size_t nb_fold, const size_t nb_class,
                          const std::vector<size_t> &label_vect, const std::vector<float> &norm_count_vect)
{
    const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::regression::LogisticRegression<> rgc(norm_count_vect, label_vect);
        MLMetrics my_f1;
        return my_f1.Evaluate(rgc, norm_count_vect, label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>, MLMetrics>
            score_data(nb_fold_final, norm_count_vect, label_vect);
        return score_data.Evaluate();
    }
}

const double CalcSVMScore(const size_t nb_fold, const size_t nb_class,
                          const std::vector<size_t> &label_vect, const std::vector<float> &norm_count_vect)
{
    mlpack::svm::LinearSVM<> lsvm(norm_count_vect, label_vect, nb_class);
    mlpack::svm::LinearSVMFunction<> lsvm_fun(norm_count_vect, label_vect, nb_class);
    return lsvm_fun.Evaluate(lsvm.Parameters());
}

#endif //KAMRAT_UTILS_FEATUREEVAL_HPP