#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

#include <cmath>
// #include <boost/math/distributions/students_t.hpp>
// #include <mlpack/core/cv/k_fold_cv.hpp>
// #include <mlpack/core/cv/metrics/f1.hpp>
// #include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
// #include <mlpack/methods/logistic_regression/logistic_regression.hpp>
// #include <mlpack/methods/linear_svm/linear_svm.hpp>
// #include <mlpack/methods/linear_svm/linear_svm_function.hpp>
// #include <boost/math/distributions/students_t.hpp>

#include "scorer.hpp" // armadillo library need be included after mlpack

/* ============================ Scoring Method ============================ *\
 * rsd          relative standard deviation                                  *
 * ttest        p-value of t-test, adjusted by Benjamini-Hochberg procedure  *
 * snr          signal-to-noise ratio                                        *
 * nbc          naive Bayes classifier                                       *
 * lrc          logistic regression classifier                               *
 * svm          support vector machine                                       *
 * user:name    use representative value in table as score                   *
\* ======================================================================== */

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
                            const double mean1, const double mean2, const double stddev1, const double stddev2)
{
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

const double CalcSNRScore(double &mean1, double &mean2, double &stddev1, double &stddev2,
                          const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect)
{
    static arma::Mat<double> norm_count_vect1 = norm_count_vect.elem(arma::find(label_vect == 1)),
                             norm_count_vect2 = norm_count_vect.elem(arma::find(label_vect == 2));
    mean1 = arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect1));
    mean2 = arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect2));
    stddev1 = arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect1));
    stddev2 = arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect2));
    if (stddev1 == 0 && stddev2 == 0)
    {
        return std::nan("");
    }
    else
    {
        return ((mean1 - mean2) / (stddev1 + stddev2));
    }
}

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

Scorer::Scorer(const ScoreMethodCode score_method_code, const std::string &score_cmd)
    : score_method_code_(score_method_code), score_cmd_(score_cmd)
{
    if (score_method_code_ == ScoreMethodCode::kNaiveBayes || score_method_code_ == ScoreMethodCode::kLogitReg)
    {
        if (score_cmd.empty())
        {
            nb_fold_ = 1; // by default, no cross-validation, use all sample for both train and test
        }
        else
        {
            nb_fold_ = std::stoul(score_cmd); // C++: typedef unsigned long size_t
        }
    }
}

const void Scorer::LoadSampleLabels(const std::vector<size_t> &label_vect)
{
    label_vect_ = arma::conv_to<arma::Row<size_t>>::from(label_vect);
    const size_t kNbClass = label_vect_.max();
    nb_smp_condi_.resize(kNbClass, 0);
    for (size_t x : label_vect)
    {
        nb_smp_condi_[x]++;
    }
    if (kNbClass != 2 && (score_method_code_ == ScoreMethodCode::kTtest ||
                          score_method_code_ == ScoreMethodCode::kSNR ||
                          score_method_code_ == ScoreMethodCode::kLog2FC ||
                          score_method_code_ == ScoreMethodCode::kLogitReg))
    {
        throw std::domain_error("T-test, signal-to-noise ratio, log2FC, and logistic regression scoring only accept binary sample condition");
    }
    if (kNbClass < 2 && (score_method_code_ == ScoreMethodCode::kNaiveBayes))
    {
        throw std::domain_error("Naive Bayes classfier only accepts condition number >= 2");
    }
}

const std::string &Scorer::GetScoreMethod() const
{
    return kScoreMethodName[score_method_code_];
}

const std::string &Scorer::GetScoreCmd() const
{
    return score_cmd_;
}

const size_t Scorer::GetNbFold() const
{
    return nb_fold_;
}

const size_t Scorer::GetNbClass() const
{
    return nb_class_;
}

const void Scorer::PrepareCountVect(const FeatureElem &feature_elem, const std::vector<double> &nf_vect, const bool to_ln)
{
    static const size_t kNbCount = label_vect_.size();
    static std::vector<double> raw_count_vect;
    norm_count_vect_.zeros(1, kNbCount);
    feature_elem.RetrieveCountVect(raw_count_vect, kNbCount);
    for (size_t i = 0; i < kNbCount; ++i)
    {
        norm_count_vect_(0, i) = (to_ln ? log(nf_vect[i] * raw_count_vect[i] + 1) : (nf_vect[i] * raw_count_vect[i]));
    }
}

const void Scorer::CalcFeatureStats(FeatureElem &feature_elem, const size_t nb_class) const
{
    feature_elem.ReserveCondiStats(nb_class);
    for (size_t i_class(1); i_class <= nb_class; ++i_class)
    {
        arma::Mat<double> &&norm_count_vect_i = norm_count_vect_.elem(arma::find(label_vect_ == i_class));
        feature_elem.AddCondiStats(arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect_i)),
                                   arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect_i)));
    }
}

const double Scorer::EvaluateFeature(FeatureElem &feature_elem, const std::vector<double> nf_vect, bool to_ln) const
{

    switch (score_method_code_)
    {
    case ScoreMethodCode::kRelatSD:
        double all_mean, all_stddev;
        feature_elem.SetScore(CalcRelatSDScore(all_mean, all_stddev, norm_count_vect));
        feature_elem.ReserveCondiStats(1);
        feature_elem.AddCondiStats(all_mean, all_stddev);
        break;
    case ScoreMethodCode::kTtest:
        double mean1, mean2, stddev1, stddev2;
        feature_elem.SetScore(CalcTtestScore(mean1, mean2, stddev1, stddev2, label_vect_, norm_count_vect));
        feature_elem.ReserveCondiStats(2);
        feature_elem.AddCondiStats(mean1, stddev1);
        feature_elem.AddCondiStats(mean2, stddev2);
        break;
    case ScoreMethodCode::kSNR:
        double mean1, mean2, stddev1, stddev2;
        feature_elem.SetScore(CalcSNRScore(mean1, mean2, stddev1, stddev2, label_vect_, norm_count_vect));
        feature_elem.ReserveCondiStats(2);
        feature_elem.AddCondiStats(mean1, stddev1);
        feature_elem.AddCondiStats(mean2, stddev2);
        break;
    }
}

#endif //KMEREVALUATE_EVALMETHODS_H
