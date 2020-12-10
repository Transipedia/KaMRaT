#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

#include <cmath>

#include "mlpack/core/cv/k_fold_cv.hpp"
#include "mlpack/core/cv/metrics/f1.hpp"
#include "mlpack/core/cv/metrics/accuracy.hpp"
#include "mlpack/methods/naive_bayes/naive_bayes_classifier.hpp"
#include "mlpack/methods/logistic_regression/logistic_regression.hpp"
#include "mlpack/methods/linear_svm/linear_svm.hpp"
#include "mlpack/methods/linear_svm/linear_svm_function.hpp"
#include "boost/math/distributions/students_t.hpp"
#include "scorer.hpp" // armadillo library need be included after mlpack

#define MET_SD "sd"
#define MET_RSD "relat.sd"
#define MET_TTEST "t-test"
#define MET_ES "effect.size"
#define MET_LFC "log2fc"
#define MET_NBF1 "naivebayes.f1"
#define MET_RGF1 "regression.f1"
#define MET_SVM "svm.hingeloss"
#define MET_USER "user"

inline float calc_stat(const arma::mat &sample_counts, const std::string &stat_name)
{
    if (stat_name == "mean")
    {
        return (arma::mean(arma::conv_to<arma::vec>::from(sample_counts)));
    }
    else if (stat_name == "median")
    {
        return (arma::median(arma::conv_to<arma::vec>::from(sample_counts)));
    }
    else if (stat_name == "sd")
    {
        return (arma::stddev(arma::conv_to<arma::vec>::from(sample_counts), 0));
    }
    else
    {
        throw std::domain_error("unknown stats name");
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

Scorer::Scorer(const std::string &score_method, const std::string &score_cmd, const std::string &sort_mode, const size_t nb_fold,
               const bool to_ln, const bool to_standardize)
    : score_method_(score_method),
      sort_mode_(sort_mode),
      score_cmd_(score_cmd),
      nb_fold_(nb_fold),
      nb_class_(0),
      to_ln_(to_ln),
      to_standardize_(to_standardize)
{
}

const std::string &Scorer::GetScoreMethod() const
{
    return score_method_;
}

const std::string &Scorer::GetSortMode() const
{
    return sort_mode_;
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

void Scorer::LoadSampleLabel(const TabHeader &tab_header)
{
    nb_class_ = tab_header.GetNbCondition();
    if (nb_class_ != 2 && (score_method_ == MET_TTEST || score_method_ == MET_ES || score_method_ == MET_LFC))
    {
        throw std::domain_error("T-test, effect size, and log2FC scoring only accept binary sample condition");
    }
    if (nb_class_ < 2 && (score_method_ == MET_NBF1 || score_method_ == MET_RGF1))
    {
        throw std::domain_error("Naive Bayes or regression classfiers only accept condition number >= 2");
    }
    std::vector<size_t> label_vect;
    tab_header.ParseSmpLabels(label_vect);
    sample_labels_ = arma::conv_to<arma::Row<size_t>>::from(label_vect); // each column is an observation
}

void Scorer::LoadSampleCount(const std::vector<float> &count_vect, const std::vector<double> &nf_vect)
{
    condi_sample_counts_.clear();
    size_t nb_sample = count_vect.size();
    sample_counts_.zeros(1, nb_sample);
    for (size_t i = 0; i < nb_sample; ++i)
    {
        sample_counts_(0, i) = nf_vect[i] * (count_vect[i] + 1);
    }
    sample_counts_.print("Normalized sample counts: ");
    for (size_t i(0); i < nb_class_; ++i)
    {
        condi_sample_counts_.emplace_back(sample_counts_.elem(arma::find(sample_labels_ == i)));
    }
}

void Scorer::ClearSampleCount()
{
    sample_counts_.clear();
    for (size_t i(0); i < condi_sample_counts_.size(); ++i)
    {
        condi_sample_counts_[i].clear();
    }
    condi_sample_counts_.clear();
}

const void Scorer::CalcCondiMeans(std::vector<float> &condi_means) const
{
    for (size_t i(0); i < nb_class_; ++i)
    {
        condi_means.emplace_back(calc_stat(condi_sample_counts_[i], "mean"));
    }
}

const void Scorer::TransformCounts() // first log then standardize
{
    size_t nb_sample = sample_counts_.size();
    if (to_ln_)
    {
        for (size_t i = 0; i < nb_sample; ++i)
        {
            sample_counts_(0, i) = log(sample_counts_(0, i)); // an offset 1 had been added on the raw count
        }
    }
    if (to_standardize_)
    {
        float mean = calc_stat(sample_counts_, "mean"), sd = calc_stat(sample_counts_, "sd");
        if (sd == 0) // in case of a constant count vector
        {
            std::cerr << "[warning]:    a constant count vector appears, turing it to all zeros." << std::endl;
            sd = 1;
        }
        for (size_t i = 0; i < nb_sample; ++i)
        {
            sample_counts_(0, i) = (sample_counts_(0, i) - mean) / sd;
        }
    }
}

const float Scorer::EvaluateScore() const // no implement, but virtual method must be defined
{
    return 0;
}

// =====> Standard Deviation Scoring <===== //
SDScorer::SDScorer(const std::string &sort_mode)
    : Scorer(MET_SD, "", (sort_mode.empty() ? "dec" : sort_mode), 1, to_ln_, to_standardize_)
{
}

const float SDScorer::EvaluateScore() const
{
    float sd = calc_stat(sample_counts_, "sd");
    return sd;
}

// =====> Relative Standard Deviation Scoring <===== //
RelatSDScorer::RelatSDScorer(const std::string &sort_mode)
    : Scorer(MET_RSD, "", (sort_mode.empty() ? "dec" : sort_mode), 1, to_ln_, to_standardize_)
{
}

const float RelatSDScorer::EvaluateScore() const
{
    float mean = calc_stat(sample_counts_, "mean"), sd = calc_stat(sample_counts_, "sd"), score = (mean <= 1 ? sd : sd / mean);
    return ((std::isnan(score) || std::isinf(score)) ? 0.0 : score);
}

// =====> T-test Scoring <===== //
TtestScorer::TtestScorer(const std::string &sort_mode)
    : Scorer(MET_TTEST, "", (sort_mode.empty() ? "inc" : sort_mode), 1, to_ln_, to_standardize_)
{
}

const float TtestScorer::EvaluateScore() const
{
    size_t cond1_num = condi_sample_counts_[0].size(), cond2_num = condi_sample_counts_[1].size();
    float cond1_mean = calc_stat(condi_sample_counts_[0], "mean"),
          cond2_mean = calc_stat(condi_sample_counts_[1], "mean"),
          cond1_sd = calc_stat(condi_sample_counts_[0], "sd"),
          cond2_sd = calc_stat(condi_sample_counts_[1], "sd"), pvalue;
    if (cond1_sd == 0 && cond2_sd == 0)
    {
        pvalue = 1;
    }
    else
    {
        float t1 = cond1_sd * cond1_sd / cond1_num,
              t2 = cond2_sd * cond2_sd / cond2_num,
              df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (cond1_num - 1) + t2 * t2 / (cond2_num - 1)),
              t_stat = (cond1_mean - cond2_mean) / sqrt(t1 + t2);
        try
        {
            boost::math::students_t dist(df);
            pvalue = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
        }
        catch (const std::exception &)
        {
            sample_counts_.print("Sample counts:");
            sample_labels_.print("Sample labels:");
            std::cout << cond1_num << "\t" << cond1_mean << "\t" << cond1_sd << "\t" << t1 << "\t" << t2 << "\t" << df << "\t" << t_stat << std::endl;
            throw std::domain_error("");
        }
    }
    return pvalue;
}

// =====> Effect Size Scorer <===== //
EffectSizeScorer::EffectSizeScorer(const std::string &sort_mode)
    : Scorer(MET_ES, "", (sort_mode.empty() ? "dec:abs" : sort_mode), 1, to_ln_, to_standardize_)
{
}

const float EffectSizeScorer::EvaluateScore() const
{
    float cond1_mean = calc_stat(condi_sample_counts_[0], "mean"),
          cond2_mean = calc_stat(condi_sample_counts_[1], "mean"),
          cond1_sd = calc_stat(condi_sample_counts_[0], "sd"),
          cond2_sd = calc_stat(condi_sample_counts_[1], "sd"),
          score = (cond2_mean - cond1_mean) / (cond1_sd + cond2_sd);
    return ((std::isnan(score) || std::isinf(score)) ? 0.0 : score);
}

// =====> Log2FC Scorer <===== //
LFCScorer::LFCScorer(const std::string &score_cmd, const std::string &sort_mode)
    : Scorer(MET_LFC, score_cmd, (sort_mode.empty() ? "dec:abs" : sort_mode), 1, to_ln_, to_standardize_)
{
}

const float LFCScorer::EvaluateScore() const
{
    float cond1_val = calc_stat(condi_sample_counts_[0], score_cmd_),
          cond2_val = calc_stat(condi_sample_counts_[1], score_cmd_),
          score = log2(cond2_val / cond1_val);
    return ((std::isnan(score) || std::isinf(score)) ? 0.0 : score);
}

// =====> Naive Bayes Scorer <===== //
NaiveBayesScorer::NaiveBayesScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer(MET_NBF1, "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold, to_ln_, to_standardize_)
{
}

// inline void CompareScore(const arma::mat &sample_counts, const arma::Row<size_t> &sample_labels, const size_t nb_class, const size_t nb_fold)
// {
//     // Binary F1 all samples
//     {
//         mlpack::naive_bayes::NaiveBayesClassifier<> nbc(sample_counts, sample_labels, nb_class);
//         arma::Row<size_t> pred_class(sample_labels.size());
//         nbc.Classify(sample_counts, pred_class);
//         mlpack::cv::F1<mlpack::cv::Binary> F1;
//         std::cout << "\t" << F1.Evaluate(nbc, sample_counts, sample_labels);
//     }
//     // Micro F1 all samples
//     {
//         mlpack::naive_bayes::NaiveBayesClassifier<> nbc(sample_counts, sample_labels, nb_class);
//         arma::Row<size_t> pred_class(sample_labels.size());
//         nbc.Classify(sample_counts, pred_class);
//         mlpack::cv::F1<mlpack::cv::Micro> F1;
//         std::cout << "\t" << F1.Evaluate(nbc, sample_counts, sample_labels);
//     }
//     // Accuracy all samples
//     {
//         mlpack::naive_bayes::NaiveBayesClassifier<> nbc(sample_counts, sample_labels, nb_class);
//         arma::Row<size_t> pred_class(sample_labels.size());
//         nbc.Classify(sample_counts, pred_class);
//         mlpack::cv::Accuracy accuracy;
//         std::cout << "\t" << accuracy.Evaluate(nbc, sample_counts, sample_labels);
//     }
//     // Binary F1 cross-validation
//     {
//         mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Binary>>
//             score_data(nb_fold, sample_counts, sample_labels, nb_class);
//         std::cout << "\t" << score_data.Evaluate();
//     }
//     // Micro F1 cross-validation
//     {
//         mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Micro>>
//             score_data(nb_fold, sample_counts, sample_labels, nb_class);
//         std::cout << "\t" << score_data.Evaluate();
//     }
//     // Accuracy cross-validation
//     {
//         mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::Accuracy>
//             score_data(nb_fold, sample_counts, sample_labels, nb_class);
//         std::cout << "\t" << score_data.Evaluate();
//     }
//     // Binary F1 cross-validation with -nan replaced by 0 in folds
//     {
//         mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
//             score_data(nb_fold, sample_counts, sample_labels, nb_class);
//         std::cout << "\t" << score_data.Evaluate() << std::endl;
//     }
// }

const float NaiveBayesScorer::EvaluateScore() const
{
    // CompareScore(sample_counts_, sample_labels_, nb_class_, nb_fold_);
    const size_t nb_fold_final = (nb_fold_ == 0 ? sample_labels_.size() : nb_fold_);
    float score;
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(sample_counts_, sample_labels_, nb_class_);
        MLMetrics my_f1;
        score = my_f1.Evaluate(nbc, sample_counts_, sample_labels_);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
            score_data(nb_fold_final, sample_counts_, sample_labels_, nb_class_);
        score = score_data.Evaluate();
    }
    if (std::isnan(score) || std::isinf(score))
    {
        std::cout << score << std::endl;
        sample_counts_.print("Sample counts: ");
        sample_labels_.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in naive Bayes model");
    }
    return score;
}

// =====> Regression Scorer <===== //
RegressionScorer::RegressionScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer(MET_RGF1, "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold, to_ln_, to_standardize_)
{
}

const float RegressionScorer::EvaluateScore() const
{
    if (nb_class_ != 2)
    {
        throw std::domain_error("Regression classifier accepts only binary condition for now");
    }
    const size_t nb_fold_final = (nb_fold_ == 0 ? sample_labels_.size() : nb_fold_);
    float score;
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::regression::LogisticRegression<> rgc(sample_counts_, sample_labels_);
        MLMetrics my_f1;
        score = my_f1.Evaluate(rgc, sample_counts_, sample_labels_);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>, MLMetrics>
            score_data(nb_fold_final, sample_counts_, sample_labels_);
        score = score_data.Evaluate();
    }
    if (std::isnan(score) || std::isinf(score))
    {
        std::cout << score << std::endl;
        sample_counts_.print("Sample counts: ");
        sample_labels_.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in regression model");
    }
    return score;
}

// =====> SVM Scorer <===== //

SVMScorer::SVMScorer(const std::string &sort_mode)
    : Scorer(MET_SVM, "", (sort_mode.empty() ? "inc" : sort_mode), 0, to_ln_, to_standardize_)
{
}

const float SVMScorer::EvaluateScore() const
{
    mlpack::svm::LinearSVM<> lsvm(sample_counts_, sample_labels_, nb_class_);
    // MLMetrics my_f1;
    // float score = my_f1.Evaluate(lsvm, sample_counts, sample_labels_);
    mlpack::svm::LinearSVMFunction<> lsvm_fun(sample_counts_, sample_labels_, nb_class_);
    float score = lsvm_fun.Evaluate(lsvm.Parameters());
    // lsvm.Parameters().print("Parameters:");
    // sample_counts.print("Sample counts:");
    // sample_labels_.print("Real labels:");
    // arma::Row<size_t> pred_class(sample_labels_.size());
    // lsvm.Classify(sample_counts, pred_class);
    // pred_class.print("Predicted labels:");
    return score;
}

// =====> User Scorer <===== //
UserScorer::UserScorer(const std::string &sort_mode)
    : Scorer(MET_USER, "", (sort_mode.empty() ? "dec" : sort_mode), 0, to_ln_, to_standardize_)
{
}

#endif //KMEREVALUATE_EVALMETHODS_H
