#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

#include <iostream>
#include "scorer.hpp"

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
    float sd = arma::stddev(arma::conv_to<arma::vec>::from(arma_sample_counts), 0);
    return sd
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

const size_t InferNbFold(const std::string &score_method, const std::string &score_cmd)
{
    if (score_method == "nb" || score_method == "sr")
    {
        return (score_cmd.empty() ? 2 : std::stoi(score_cmd));
    }
    else
    {
        return 0;
    }
}

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
    if (nb_class != 2 && (score_method_ == "ttest" || score_method_ == "es" || score_method_ == "lfc"))
    {
        throw std::domain_error("T-test, effect size, and log2FC scoring only accept binary sample condition");
    }
    nb_class_ = nb_class;
    sample_labels_ = arma::conv_to<arma::Row<size_t>>::from(label_vect); // each column is an observation
    sample_labels_.print("Sample labels: ");
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

// =====> Standard Deviation Scoring <===== //
SDScorer::SDScorer(const std::string &sort_mode)
    : Scorer("sd", "", (sort_mode.empty() ? "dec" : sort_mode), 1)
{
}

const float SDScorer::CalcScore(const std::vector<float> &sample_counts) const override
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    float sd = calc_sd(arma_sample_counts);
    return sd;
}

// =====> Relative Standard Deviation Scoring <===== //
RelatSDScorer::RelatSdScorer(const std::string &sort_mode)
    : Scorer("relat.sd", "", (sort_mode.empty() ? "dec" : sort_mode), 1)
{
}

const float RelatSDScorer::CalcScore(const std::vector<float> &sample_counts) const override
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
    : Scorer("t-test", "", (sort_mode.empty() ? "inc" : sort_mode), 1)
{
}

const float TtestScorer::CalcScore(const std::vector<float> &sample_counts) const override
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
EffectSizeScorer::EffectSizeScorer(const std::string sort_mode)
    : Scorer("effect size", "", (sort_mode.empty() ? "dec:abs" : sort_mode), 1)
{
}

const float EffectSizeScorer::CalcScore(const std::vector<float> &sample_counts) const override
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
    : Scorer("log2fc", score_cmd, (sort_mode.empty() ? "dec:abs" : sort_mode), 1)
{
}

const float LFCScorer::CalcScore(const std::vector<float> &sample_counts) const override
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
        throw std::domain_error("unknown log2FC command: " + score_cmd);
    }
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

// =====> Naive Bayes Scorer <===== //
NaiveBayesScorer::NaiveBayesScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer("naive bayes f1 score", "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold)
{
}

const float NaiveBayesScorer::CalcScore(const std::vector<float> &sample_counts) const override
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    float score;
    if (nb_class_ == 2)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Binary>>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        score = score_data.Evaluate();
    }
    else if (nb_class > 2)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Micro>>
            score_data(nb_fold, arma_sample_counts, sample_labels, nb_class);
        score = score_data.Evaluate();
    }
    else
    {
        throw std::invalid_argument("nb_class should be at least 2 for naive bayes classifier (nb_class = " + std::to_string(nb_class_) + ")");
    }
    if (isnan(score) || isinf(score))
    {
        arma_sample_counts.print("Sample counts: ");
        arma_sample_counts.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in naive Bayes model");
    }
    return score;
}

// =====> Regression Scorer <===== //
RegressionScorer::RegressionScorer(const std::string &sort_mode, const size_t nb_fold)
    : Scorer("regression f1 score", "", (sort_mode.empty() ? "dec" : sort_mode), nb_fold)
{
}

const float RegressionScorer::CalcScore(const std::vector<float> &sample_counts) const override
{
    arma::mat arma_sample_counts;
    Conv2Arma(arma_sample_counts, sample_counts);
    float score;
    if (nb_class_ == 2)
    {
        mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression, mlpack::cv::F1<mlpack::cv::Binary>>
            score_data(nb_fold, sample_counts, sample_labels, nb_class);
        score = score_data.Evaluate();
    }
    else if (nb_class_ > 2)
    {
        mlpack::cv::KFoldCV<mlpack::regression::SoftmaxRegression, mlpack::cv::F1<mlpack::cv::Micro>>
            score_data(nb_fold, sample_counts, sample_labels, nb_class);
        score = score_data.Evaluate();
    }
    else
    {
        throw std::invalid_argument("nb_class should be at least 2 for regression classifier (nb_class = " + std::to_string(nb_class_) + ")");
    }
    if (isnan(score) || isinf(score))
    {
        sample_counts.print("Sample counts: ");
        sample_labels.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in regression model");
    }
    return score;
}

// =====> User Scorer <===== //
UserScorer::UserScorer(const std::string &sort_mode)
    : Scorer("user", "", (sort_mode.empty() ? "dec" : sort_mode), 0)
{
}

// const float UserScorer::CalcScore(const std::vector<float> &sample_counts) const override // no implement, 
// {
//     return 0;
// }

#endif //KMEREVALUATE_EVALMETHODS_H
