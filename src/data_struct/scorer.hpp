#ifndef KAMRAT_DATASTRUCT_SCORER_HPP
#define KAMRAT_DATASTRUCT_SCORER_HPP

#include <string>
#include <vector>
#include <cmath>

#include "armadillo"
#include "tab_header.hpp"

class Scorer
{
public:
    Scorer(const std::string &score_method, const std::string &score_cmd, const std::string &sort_mode, size_t nb_fold);
    void LoadSampleLabel(const TabHeader &tab_header);
    const std::string &GetScoreMethod() const;
    const std::string &GetSortMode() const;
    const std::string &GetScoreCmd() const;
    const size_t GetNbFold() const;
    const size_t GetNbClass() const;
    virtual const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const;

protected:
    const std::string score_method_;  // nb, lr, sd, rsd, ttest, es, lfc, user
    const std::string sort_mode_;     // dec, dec::abs, inc, inc::abs
    const std::string score_cmd_;     // to nb_fold_ if nb or lr; else for mean or median in lfc or for score colname in user
    const size_t nb_fold_;            // for naive Bayes and logistic regression
    size_t nb_class_;                 // number of conditions
    arma::Row<size_t> sample_labels_; // sample labels
};

class SDScorer : public Scorer
{
public:
    SDScorer(const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class RelatSDScorer : public Scorer
{
public:
    RelatSDScorer(const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class TtestScorer : public Scorer
{
public:
    TtestScorer(const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class EffectSizeScorer : public Scorer
{
public:
    EffectSizeScorer(const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class LFCScorer : public Scorer
{
public:
    LFCScorer(const std::string &score_cmd, const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class NaiveBayesScorer : public Scorer
{
public:
    NaiveBayesScorer(const std::string &sort_mode, size_t nb_fold);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class RegressionScorer : public Scorer
{
public:
    RegressionScorer(const std::string &sort_mode, size_t nb_fold);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class SVMScorer : public Scorer
{
public:
    SVMScorer(const std::string &sort_mode);
    const float CalcScore(const std::vector<float> &sample_counts, const bool to_standardize) const override;
};

class UserScorer : public Scorer
{
public:
    UserScorer(const std::string &sort_mode);
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
