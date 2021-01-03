#ifndef KAMRAT_DATASTRUCT_SCORER_HPP
#define KAMRAT_DATASTRUCT_SCORER_HPP

#include <string>
#include <vector>
#include <armadillo>

#include "feature_elem.hpp"

enum ScoreMethodCode
{
    kRelatSD = 0,
    kTtest,
    kSNR,
    kLog2FC,
    kNaiveBayes,
    kLogitReg,
    kSVM,
    kUser
};
const std::vector<std::string> kScoreMethodName{"relat.sd", "t-test", "SNR", "log2fc", "naive.bayes", "logit.reg", "SVM", "user"};

class Scorer
{
public:
    Scorer(const ScoreMethodCode score_method_code, const std::string &score_cmd);

    const void LoadSampleLabels(const std::vector<size_t> &label_vect);

    const std::string &GetScoreMethod() const;
    const std::string &GetScoreCmd() const;
    const size_t GetNbFold() const;
    const size_t GetNbClass() const;

    const void PrepareCountVect(const FeatureElem &feature_elem, const std::vector<double> &nf_vect, bool to_ln);
    const void CalcFeatureStats(FeatureElem &feature_elem, size_t nb_class) const;
    const double EvaluateFeature(FeatureElem &feature_elem) const;

private:
    const ScoreMethodCode score_method_code_; // scoring method code
    const std::string score_cmd_;             // fold number if naive Bayes or logit regression; mean or median if lfc; or score colname if user
    size_t nb_fold_;                          // number of fold for naive Bayes or logit regression classifiers
    size_t nb_class_;                         // number of conditions for naive Bayes or SVM
    std::vector<size_t> nb_smp_condi_;        // number of sample in each condition
    arma::Row<size_t> label_vect_;            // sample label vector
    arma::Mat<double> norm_count_vect_;       // normalized (and transformed) sample count vector
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
