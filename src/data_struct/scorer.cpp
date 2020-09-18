#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

#include <iostream>
#include "scorer.hpp"

#define MET_SD "sd"
#define MET_RSD "relat.sd"
#define MET_TTEST "t-test"
#define MET_ES "effect.size"
#define MET_LFC "log2fc"
#define MET_NBF1 "naivebayes.f1"
#define MET_RGF1 "regression.f1"
#define MET_USER "user"

inline float calc_mean(const arma::mat &sample_counts)
{
    float mean = arma::mean(arma::conv_to<arma::vec>::from(sample_counts));
    return mean;
}

inline float calc_median(const arma::mat &sample_counts)
{
    float median = arma::median(arma::conv_to<arma::vec>::from(sample_counts));
    return median;
}

inline float calc_sd(const arma::mat &sample_counts)
{
    float sd = arma::stddev(arma::conv_to<arma::vec>::from(sample_counts), 0);
    return sd;
}

inline void Conv2Arma(arma::mat &arma_count_vect, const std::vector<float> &count_vect)
{
    size_t nb_sample = count_vect.size();
    arma_count_vect.zeros(1, nb_sample);
    for (int i = 0; i < nb_sample; ++i)
    {
        arma_count_vect(0, i) = count_vect[i];
    }
}

class MLMetrics
{
public:
    template <typename MLAlgorithm, typename DataType>
    static double Evaluate(MLAlgorithm &model, const DataType &data, const arma::Row<size_t> &labels)
    {
        arma::Row<size_t> pred_class(labels.size());
        model.Classify(data, pred_class);
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

Scorer::Scorer(const std::string &score_method, const std::string &score_cmd, const std::string &sort_mode, const size_t nb_fold)
    : score_method_(score_method),
      sort_mode_(sort_mode),
      score_cmd_(score_cmd),
      nb_fold_(nb_fold),
      nb_class_(0)
{
}

void Scorer::LoadSampleLabel(const std::vector<size_t> &label_vect, const size_t nb_class)
{
    if (nb_class != 2 && (score_method_ == MET_TTEST || score_method_ == MET_ES || score_method_ == MET_LFC))
    {
        throw std::domain_error("T-test, effect size, and log2FC scoring only accept binary sample condition");
    }
    if (nb_class < 2 && (score_method_ == MET_NBF1 || score_method_ == MET_RGF1))
    {
        throw std::domain_error("Naive Bayes or regression classfiers only accept condition number >= 2");
    }
    nb_class_ = nb_class;
    sample_labels_ = arma::conv_to<arma::Row<size_t>>::from(label_vect); // each column is an observation
    // sample_labels_.print("Sample labels: ");
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

const float Scorer::CalcScore(const std::vector<float> &sample_counts) const // no implement, but virtual method must be defined
{
    return 0;
}

// =====> Standard Deviation Scoring <===== //
SDScorer::SDScorer(const std::string &sort_mode)
    : Scorer(MET_SD, "", (sort_mode.empty() ? "dec" : sort_mode), 1)
{
}

const float SDScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    float sd = calc_sd(arma_sample_counts);
    return sd;
}

// =====> Relative Standard Deviation Scoring <===== //
RelatSDScorer::RelatSDScorer(const std::string &sort_mode)
    : Scorer(MET_RSD, "", (sort_mode.empty() ? "dec" : sort_mode), 1)
{
}

const float RelatSDScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    float sd = calc_sd(arma_sample_counts),
          mean = calc_mean(arma_sample_counts),
          score = (mean <= 1 ? sd : sd / mean);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

// =====> T-test Scoring <===== //
TtestScorer::TtestScorer(const std::string &sort_mode)
    : Scorer(MET_TTEST, "", (sort_mode.empty() ? "inc" : sort_mode), 1)
{
}

const float TtestScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    arma::mat cond1_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 0)),
              cond2_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 1));
    for (size_t i = 0; i < cond1_counts.size(); ++i)
    {
        cond1_counts(i, 0) = log(cond1_counts(i, 0) + 1);
    }
    for (size_t i = 0; i < cond2_counts.size(); ++i)
    {
        cond2_counts(i, 0) = log(cond2_counts(i, 0) + 1);
    }
    size_t cond1_num = cond1_counts.size(), cond2_num = cond2_counts.size();
    float cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts),
          cond1_sd = calc_sd(cond1_counts), cond2_sd = calc_sd(cond2_counts),
          pvalue;
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
        boost::math::students_t dist(df);
        pvalue = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
    }
    return pvalue;
}

// =====> Effect Size Scorer <===== //
EffectSizeScorer::EffectSizeScorer(const std::string &sort_mode)
    : Scorer(MET_ES, "", (sort_mode.empty() ? "dec:abs" : sort_mode), 1)
{
}

const float EffectSizeScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    arma::mat cond1_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 0)),
              cond2_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 1));
    float cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts),
          cond1_sd = calc_sd(cond1_counts), cond2_sd = calc_sd(cond2_counts),
          score = (cond2_mean - cond1_mean) / (cond1_sd + cond2_sd);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

// =====> Log2FC Scorer <===== //
LFCScorer::LFCScorer(const std::string &score_cmd, const std::string &sort_mode)
    : Scorer(MET_LFC, score_cmd, (sort_mode.empty() ? "dec:abs" : sort_mode), 1)
{
}

const float LFCScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    arma::mat cond1_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 0)),
              cond2_counts = arma_sample_counts.elem(arma::find(sample_labels_ == 1));
    float score;
    if (score_cmd_ == "mean")
    {
        float cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts);
        score = log2(cond2_mean / cond1_mean);
    }
    else if (score_cmd_ == "median")
    {
        float cond1_median = calc_median(cond1_counts), cond2_median = calc_median(cond2_counts);
        score = log2(cond2_median / cond1_median);
    }
    else
    {
        throw std::domain_error("unknown log2FC command: " + score_cmd_);
    }
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

// =====> Naive Bayes Scorer <===== //
NaiveBayesScorer::NaiveBayesScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer(MET_NBF1, "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold)
{
}

inline void CompareScore(const arma::mat &arma_sample_counts, const arma::Row<size_t> &sample_labels, const size_t nb_class, const size_t nb_fold)
{
    // Binary F1 all samples
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(arma_sample_counts, sample_labels, nb_class);
        arma::Row<size_t> pred_class(sample_labels.size());
        nbc.Classify(arma_sample_counts, pred_class);
        mlpack::cv::F1<mlpack::cv::Binary> F1;
        std::cout << "\t" << F1.Evaluate(nbc, arma_sample_counts, sample_labels);
    }
    // Micro F1 all samples
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(arma_sample_counts, sample_labels, nb_class);
        arma::Row<size_t> pred_class(sample_labels.size());
        nbc.Classify(arma_sample_counts, pred_class);
        mlpack::cv::F1<mlpack::cv::Micro> F1;
        std::cout << "\t" << F1.Evaluate(nbc, arma_sample_counts, sample_labels);
    }
    // Accuracy all samples
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(arma_sample_counts, sample_labels, nb_class);
        arma::Row<size_t> pred_class(sample_labels.size());
        nbc.Classify(arma_sample_counts, pred_class);
        mlpack::cv::Accuracy accuracy;
        std::cout << "\t" << accuracy.Evaluate(nbc, arma_sample_counts, sample_labels);
    }
    // Binary F1 cross-validation
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Binary>>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        std::cout << "\t" << score_data.Evaluate();
    }
    // Micro F1 cross-validation
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Micro>>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        std::cout << "\t" << score_data.Evaluate();
    }
    // Accuracy cross-validation
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::Accuracy>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        std::cout << "\t" << score_data.Evaluate();
    }
    // Binary F1 cross-validation with -nan replaced by 0 in folds
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        std::cout << "\t" << score_data.Evaluate() << std::endl;
    }
}

const float NaiveBayesScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    // CompareScore(arma_sample_counts, sample_labels_, nb_class_, nb_fold_);
    const size_t nb_fold_final = (nb_fold_ == 0 ? sample_labels_.size() : nb_fold_);
    float score;
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(arma_sample_counts, sample_labels_, nb_class_);
        arma::Row<size_t> pred_class(sample_counts.size());
        nbc.Classify(arma_sample_counts, pred_class);
        MLMetrics F1;
        score = F1.Evaluate(nbc, arma_sample_counts, sample_labels_);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
            score_data(nb_fold_final, arma_sample_counts, sample_labels_, nb_class_);
        score = score_data.Evaluate();
    }
    if (isnan(score) || isinf(score))
    {
        std::cout << score << std::endl;
        arma_sample_counts.print("Sample counts: ");
        sample_labels_.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in naive Bayes model");
    }
    return score;
}

// =====> Regression Scorer <===== //
RegressionScorer::RegressionScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer(MET_RGF1, "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold)
{
}

const float RegressionScorer::CalcScore(const std::vector<float> &sample_counts) const
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    const size_t nb_fold_final = (nb_fold_ == 0 ? sample_labels_.size() : nb_fold_);
    float score;
    if (nb_fold_final == 1)
    {
        // mlpack::regression::SoftmaxRegression rgc(arma_sample_counts, sample_labels_, nb_class_);
        // arma::Row<size_t> pred_class(sample_counts.size());
        // rgc.Classify(arma_sample_counts, pred_class);
        // arma::mat confusion_mat;
        // mlpack::data::ConfusionMatrix(pred_class, sample_labels_, confusion_mat, nb_class_);
        // pred_class.print("Pred class:");
        // sample_labels_.print("Sample labels:");
        // confusion_mat.print("Confusion matrix:");
        // MLMetrics F1;
        // score = F1.Evaluate(rgc, arma_sample_counts, sample_labels_);
    }
    else
    {
        // mlpack::cv::KFoldCV<mlpack::regression::SoftmaxRegression, MLMetrics>
        //     score_data(nb_fold_, arma_sample_counts, sample_labels_, nb_class_);
        // score = score_data.Evaluate();
    }
    throw std::domain_error("Softmax regression is not applicable for now");
    if (isnan(score) || isinf(score))
    {
        std::cout << score << std::endl;
        arma_sample_counts.print("Sample counts: ");
        sample_labels_.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in regression model");
    }
    return score;
}

// =====> User Scorer <===== //
UserScorer::UserScorer(const std::string &sort_mode)
    : Scorer(MET_USER, "", (sort_mode.empty() ? "dec" : sort_mode), 0)
{
}

#endif //KMEREVALUATE_EVALMETHODS_H
