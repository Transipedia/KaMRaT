#ifndef KMEREVALUATE_EVALMETHODS_H
#define KMEREVALUATE_EVALMETHODS_H

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

inline float sd_scoring(const arma::mat &sample_counts)
{
    float sd = arma::stddev(arma::conv_to<arma::vec>::from(sample_counts), 0);
    return sd;
}

inline float relatsd_scoring(const arma::mat &sample_counts)
{
    float sd = sd_scoring(sample_counts),
          mean = calc_mean(sample_counts),
          score = (mean <= 1 ? sd : sd / mean);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

inline double naivebayes_scoring(const arma::mat &sample_counts,
                                 const arma::Row<size_t> &sample_labels,
                                 const size_t nb_fold,
                                 const size_t nb_class)
{
    mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::F1<mlpack::cv::Micro>>
        score_data(nb_fold, sample_counts, sample_labels, nb_class);
    double score = score_data.Evaluate();
    if (isnan(score) || isinf(score))
    {
        sample_counts.print("Sample counts: ");
        sample_labels.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in naive Bayes");
    }
    return score;
}

inline double softmaxreg_scoring(const arma::mat &sample_counts,
                                 const arma::Row<size_t> &sample_labels,
                                 const size_t nb_fold,
                                 const size_t nb_class)
{
    mlpack::cv::KFoldCV<mlpack::regression::SoftmaxRegression, mlpack::cv::F1<mlpack::cv::Micro>>
        score_data(nb_fold, sample_counts, sample_labels, nb_class);
    double score = score_data.Evaluate();
    if (isnan(score) || isinf(score))
    {
        sample_counts.print("Sample counts: ");
        sample_labels.print("Sample labels: ");
        throw std::domain_error("F1 score is NaN or Inf in logistic regression");
    }
    return score;
}

inline double ttest_scoring(const arma::mat &sample_counts,
                            const arma::Row<size_t> &sample_labels) // mean, std, slightly different from dekupl ttestFilter, but same with R
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0)),
              cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    for (size_t i = 0; i < cond1_counts.size(); ++i)
    {
        cond1_counts(i, 0) = log(cond1_counts(i, 0) + 1);
    }
    for (size_t i = 0; i < cond2_counts.size(); ++i)
    {
        cond2_counts(i, 0) = log(cond2_counts(i, 0) + 1);
    }
    size_t cond1_num = cond1_counts.size(), cond2_num = cond2_counts.size();
    double cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts),
           cond1_sd = sd_scoring(cond1_counts), cond2_sd = sd_scoring(cond2_counts);
    double pvalue;
    if (cond1_sd == 0 && cond2_sd == 0)
    {
        pvalue = 1;
    }
    else
    {
        double t1 = cond1_sd * cond1_sd / cond1_num,
               t2 = cond2_sd * cond2_sd / cond2_num,
               df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (cond1_num - 1) + t2 * t2 / (cond2_num - 1)),
               t_stat = (cond1_mean - cond2_mean) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        pvalue = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
    }
    return pvalue;
}

inline double effectsize_scoring(const arma::mat &sample_counts,
                                 const arma::Row<size_t> &sample_labels)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0)),
              cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts),
           cond1_sd = sd_scoring(cond1_counts), cond2_sd = sd_scoring(cond2_counts),
           score = (cond2_mean - cond1_mean) / (cond1_sd + cond2_sd);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

inline double lfcmean_scoring(const arma::mat &sample_counts,
                              const arma::Row<size_t> &sample_labels)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0)),
              cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_mean = calc_mean(cond1_counts), cond2_mean = calc_mean(cond2_counts),
           score = log2(cond2_mean / cond1_mean);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

inline double lfcmedian_scoring(const arma::mat &sample_counts,
                                const arma::Row<size_t> &sample_labels)
{
    arma::mat cond1_counts = sample_counts.elem(arma::find(sample_labels == 0)),
              cond2_counts = sample_counts.elem(arma::find(sample_labels == 1));
    double cond1_median = calc_median(cond1_counts), cond2_median = calc_median(cond2_counts),
           score = log2(cond2_median / cond1_median);
    return ((isnan(score) || isinf(score)) ? 0.0 : score);
}

inline double lfc_scoring(const arma::mat &sample_counts,
                          const arma::Row<size_t> &sample_labels,
                          const std::string &lfc_cmd)
{
    if (lfc_cmd == "mean")
    {
        return lfcmean_scoring(sample_counts, sample_labels);
    }
    else if (lfc_cmd == "median")
    {
        return lfcmedian_scoring(sample_counts, sample_labels);
    }
    else
    {
        std::domain_error("unknown log2FC command" + lfc_cmd);
    }
}

const std::string &&InferSortMode(const std::string &score_method, const std::string &sort_mode)
{
    if (score_method == "sd")
    {
        return (sort_mode.empty() ? "dec" : sort_mode);
    }
    else if (score_method == "rsd")
    {
        return (sort_mode.empty() ? "dec" : sort_mode);
    }
    else if (score_method == "nb")
    {
        return (sort_mode.empty() ? "dec" : sort_mode);
    }
    else if (score_method == "lr")
    {
        return (sort_mode.empty() ? "dec" : sort_mode);
    }
    else if (score_method == "ttest")
    {
        return (sort_mode.empty() ? "inc" : sort_mode);
    }
    else if (score_method == "es")
    {
        return (sort_mode.empty() ? "dec:abs" : sort_mode);
    }
    else if (score_method == "lfc")
    {
        return (sort_mode.empty() ? "dec:abs" : sort_mode);
    }
    else if (score_method == "user")
    {
        return (sort_mode.empty() ? "dec" : sort_mode);
    }
    else
    {
        throw std::domain_error("unknown scoring method name " + score_method);
    }
}

Scorer::Scorer(const std::string &score_method, const std::string &score_cmd, const std::string &sort_mode)
    : score_method_(score_method),
      sort_mode_(InferSortMode(score_method, sort_mode)),
      score_cmd_(score_cmd),
      nb_fold_(0),
      nb_class_(0)
{
}

void Scorer::LoadSampleLabel(const std::vector<int> &label_vect, const size_t nb_class)
{
    if (score_method_ == "nb" || score_method_ == "lr")
    {
        nb_fold_ = (score_cmd_.empty() ? 2 : std::stoi(score_cmd_));
    }
    if (nb_class != 2 && (score_method_ == "ttest" || score_method_ == "es" || score_method_ == "lfc"))
    {
        throw std::domain_error("Ttest, effect size, and log2FC scoring only accept binary sample condition");
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

const float Scorer::CalcScore(const std::vector<float> &count_vect) const
{
    size_t nb_sample = count_vect.size();
    arma::mat sample_counts(1, nb_sample); // a matrix-formed row vector, each column is an observation
    for (int i = 0; i < nb_sample; ++i)
    {
        sample_counts(0, i) = count_vect[i];
    }
    if (score_method_ == "sd")
    {
        return sd_scoring(sample_counts);
    }
    else if (score_method_ == "rsd")
    {
        return relatsd_scoring(sample_counts);
    }
    else if (score_method_ == "nb")
    {
        return naivebayes_scoring(sample_counts, sample_labels_, nb_fold_, nb_class_);
    }
    else if (score_method_ == "lr")
    {
        return softmaxreg_scoring(sample_counts, sample_labels_, nb_fold_, nb_class_);
    }
    else if (score_method_ == "ttest")
    {
        return ttest_scoring(sample_counts, sample_labels_);
    }
    else if (score_method_ == "es")
    {
        return effectsize_scoring(sample_counts, sample_labels_);
    }
    else if (score_method_ == "lfc")
    {
        return lfc_scoring(sample_counts, sample_labels_, score_cmd_);
    }
    else
    {
        throw std::domain_error("unknown scoring method");
    }
}

#endif //KMEREVALUATE_EVALMETHODS_H