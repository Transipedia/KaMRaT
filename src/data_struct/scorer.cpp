// #include <boost/math/distributions/students_t.hpp>
// #include <mlpack/core/cv/k_fold_cv.hpp>
// #include <mlpack/core/cv/metrics/f1.hpp>
// #include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
// #include <mlpack/methods/logistic_regression/logistic_regression.hpp>
// #include <mlpack/methods/linear_svm/linear_svm.hpp>
// #include <mlpack/methods/linear_svm/linear_svm_function.hpp>

#include "scorer.hpp" // armadillo library need be included after mlpack

/* ============================ Scoring Method ============================ *\
 * ttest        p-value of t-test, adjusted by Benjamini-Hochberg procedure  *
 * snr          signal-to-noise ratio                                        *
 * lr           logistic regression                                          *
 * nbc          naive Bayes classifier                                       *
 * svm          support vector machine                                       *
\* ======================================================================== */

const ScorerCode ParseScorerCode(const std::string &scorer_str)
{
    if (scorer_str == "ttest")
    {
        return ScorerCode::kTtest;
    }
    else if (scorer_str == "snr")
    {
        return ScorerCode::kSNR;
    }
    else if (scorer_str == "lr")
    {
        return ScorerCode::kLR;
    }
    else if (scorer_str == "nbc")
    {
        return ScorerCode::kNBC;
    }
    else if (scorer_str == "svm")
    {
        return ScorerCode::kSVM;
    }
}

// const double CalcRelatSDScore(const double all_mean, const double all_stddev)
// {
//     if (all_mean <= 1)
//     {
//         return all_stddev;
//     }
//     else
//     {
//         return (all_stddev / all_mean);
//     }
// }

const double CalcTtestScore(const arma::Mat<float> &arma_count_vect1, const arma::Mat<float> &arma_count_vect2)
{
    double mean1 = arma::mean(arma::mean(arma_count_vect1)), mean2 = arma::mean(arma::mean(arma_count_vect2)),
           stddev1 = arma::mean(arma::stddev(arma_count_vect1)), stddev2 = arma::mean(arma::stddev(arma_count_vect2));
    size_t nb1 = arma_count_vect1.size(), nb2 = arma_count_vect2.size();
    if (stddev1 == 0 && stddev2 == 0)
    {
        return 1;
    }
    else
    {
        double t1 = stddev1 * stddev1 / nb1,
               t2 = stddev2 * stddev2 / nb2,
               df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (nb1 - 1) + t2 * t2 / (nb2 - 1)),
               t_stat = (mean1 - mean2) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        return (2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat))));
    }
}

// const double CalcSNRScore(const double mean1, const double mean2, const double stddev1, const double stddev2)
// {
//     if (stddev1 == 0 && stddev2 == 0)
//     {
//         return std::nan("");
//     }
//     else
//     {
//         return ((mean1 - mean2) / (stddev1 + stddev2));
//     }
// }

// class MLMetrics
// {
// public:
//     template <typename MLAlgorithm, typename DataType>
//     static double Evaluate(MLAlgorithm &model, const DataType &data, const arma::Row<size_t> &labels)
//     {
//         if (labels.max() == 1) // binary conditions
//         {
//             // data.print("Data:");
//             double score = mlpack::cv::F1<mlpack::cv::Binary>().Evaluate(model, data, labels);
//             return (isnan(score) ? 0 : score);
//         }
//         else // multiple conditions
//         {
//             double score = mlpack::cv::F1<mlpack::cv::Micro>().Evaluate(model, data, labels);
//             return score;
//         }
//     }
// };

// const double CalcLRCScore(const size_t nb_fold, const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect)
// {
//     const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
//     if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
//     {
//         mlpack::regression::LogisticRegression<> rgc(norm_count_vect, label_vect);
//         MLMetrics my_f1;
//         return my_f1.Evaluate(rgc, norm_count_vect, label_vect);
//     }
//     else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
//     {
//         mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>, MLMetrics>
//             score_data(nb_fold_final, norm_count_vect, label_vect);
//         return score_data.Evaluate();
//     }
// }

// const double CalcNBCScore(const size_t nb_fold, const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect, const size_t nb_class)
// {
//     const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
//     if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
//     {
//         mlpack::naive_bayes::NaiveBayesClassifier<> nbc(norm_count_vect, label_vect, nb_class);
//         MLMetrics my_f1;
//         return my_f1.Evaluate(nbc, norm_count_vect, label_vect);
//     }
//     else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
//     {
//         mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
//             score_data(nb_fold_final, norm_count_vect, label_vect, nb_class);
//         return score_data.Evaluate();
//     }
// }

// const double CalcSVMScore(const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect, const size_t nb_class)
// {
//     mlpack::svm::LinearSVM<> lsvm(norm_count_vect, label_vect, nb_class);
//     mlpack::svm::LinearSVMFunction<> lsvm_fun(norm_count_vect, label_vect, nb_class);
//     return lsvm_fun.Evaluate(lsvm.Parameters());
// }

Scorer::Scorer(const std::string &scorer_str, const size_t nfold, const std::vector<double>)
    : scorer_code_(ParseScorerCode(scorer_str)), nfold_(nfold)
{
}

const std::string &Scorer::GetScorerName() const
{
    return kScorerNameVect[scorer_code_];
}

const void Scorer::LoadSampleLabels(const std::vector<size_t> &label_vect)
{
    arma_label_vect_ = arma::conv_to<arma::Row<size_t>>::from(label_vect);
    // arma_label_vect_.print("Label vector:");
    nclass_ = arma_label_vect_.max() + 1;
    if (nclass_ != 2 && (scorer_code_ == ScorerCode::kTtest || scorer_code_ == ScorerCode::kSNR || scorer_code_ == ScorerCode::kLR))
    {
        throw std::domain_error("scoring by T-test, signal-to-noise ratio, and logistic regression only accept binary sample condition");
    }
    if (nclass_ < 2 && (scorer_code_ == ScorerCode::kNBC || scorer_code_ == ScorerCode::kSVM))
    {
        throw std::domain_error("scoring by naive Bayes or SVM only accepts condition number >= 2");
    }
}

const double Scorer::EstimateScore(const std::vector<float> &count_vect, const bool no_norm, const bool ln_transf, const bool standardize) const
{
    static arma::Mat<float> arma_count_vect = arma::conv_to<arma::Mat<float>>::from(count_vect);
    arma_count_vect.print("Count vector before transformation: ");
    if (!no_norm)
    {
        arma_count_vect = dot(arma_count_vect, arma_nf_vect_);
    }
    if (ln_transf)
    {
        arma_count_vect = log(arma_count_vect);
    }
    if (standardize)
    {
        arma_count_vect = (arma_count_vect - arma::mean(arma::mean(arma_count_vect))) / arma::mean(arma::stddev(arma_count_vect));
    }
    arma_count_vect.print("Count vector after transformation:");

    switch (scorer_code_)
    {
    case ScorerCode::kTtest:
        return CalcTtestScore(arma_count_vect);
    case ScorerCode::kSNR:
        return CalcSNRScore(feature_elem.GetCondiMeanAt(0), feature_elem.GetCondiMeanAt(1),
                            feature_elem.GetCondiStddevAt(0), feature_elem.GetCondiStddevAt(1));
    case ScorerCode::kLR:
        return CalcLRCScore(nb_fold_, arma_label_vect_, norm_count_vect_);
    case ScorerCode::kNBC:
        return CalcNBCScore(nb_fold_, arma_label_vect_, norm_count_vect_, nclass_);
    case ScorerCode::kSVM:
        return CalcSVMScore(arma_label_vect_, norm_count_vect_, nclass_);
    }
}
