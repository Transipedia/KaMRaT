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
    kLogitReg,
    kNaiveBayes,
    kSVM,
    kUser
};
const std::vector<std::string> kScoreMethodName{"relat.sd", "t-test", "SNR", "logit.reg", "naive.bayes", "SVM", "user"};
enum SortModeCode
{
    kInc = 0,
    kIncAbs,
    kDec,
    kDecAbs
};
const std::vector<std::string> kSortModeName{"ascending order",
                                             "ascending order of absolute value",
                                             "descending order",
                                             "descending order of absolute value"};

class Scorer
{
public:
    Scorer(ScoreMethodCode);                                             // relat.sd, t-test, snr, svm
    Scorer(ScoreMethodCode score_method_code, const size_t nb_fold);     // naive bayes, logit regression
    Scorer(const std::string &rep_colname, SortModeCode sort_mode_code); // column name given by user

    const void LoadSampleLabels(const std::vector<size_t> &label_vect);

    const ScoreMethodCode GetScoreMethodCode() const;
    const std::string &GetRepColname() const;
    const SortModeCode GetSortModeCode() const;
    const size_t GetNbFold() const;
    const size_t GetNbClass() const;

    const void PrepareCountVect(const FeatureElem &feature_elem, const std::vector<double> &nf_vect,
                                bool to_ln, bool to_standardize, bool no_norm);
    const void CalcFeatureStats(FeatureElem &feature_elem);
    const double EvaluateScore(const FeatureElem &feature_elem) const;

private:
    const ScoreMethodCode score_method_code_;       // scoring method code
    const std::string rep_colname_;                 // score colname given by user (for user-defined ranking)
    const SortModeCode sort_mode_code_;             // sorting mode code
    const size_t nb_fold_;                          // number of fold (for naive Bayes or logit regression classifiers)
    size_t nb_class_;                               // number of conditions for naive Bayes or SVM
    arma::Row<size_t> label_vect_;                  // sample label vector
    arma::Mat<double> norm_count_vect_;             // normalized (and transformed) sample count vector
    std::vector<size_t> nbsmp_condi_;               // number of sample in each condition
    std::vector<double> mean_condi_, stddev_condi_; // mean and standard deviation of sample count in each condition
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
