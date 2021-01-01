#ifndef KAMRAT_DATASTRUCT_SCORER_HPP
#define KAMRAT_DATASTRUCT_SCORER_HPP

#include <string>
#include <vector>
#include <cmath>

#include <armadillo>

class Scorer
{
public:
    Scorer(const std::string &score_method);

    void LoadSampleLabel(const std::vector<size_t> &label_vect);

    const std::string &GetScoreMethod() const;
    const size_t GetNbFold() const;

    const double EvaluateFeature(std::vector<std::tuple<double, double, double>> &condi_mms,
                                 const std::vector<float> &norm_count_vect) const;

private:
    enum ScoreMethod
    {
        kRelatSD,
        kTtest,
        kSNR,
        kLog2FC,
        kNaiveBayes,
        kLogitReg,
        kSVM,
        kUser
    };
    const ScoreMethod score_method_code_; // score method, coded as number
    const std::string score_cmd_;         // to nb_fold_ if nb or lr; else for mean or median in lfc or for score colname in user
    const size_t nb_fold_;                // for naive Bayes and logistic regression
    size_t nb_class_;                     // number of conditions
    arma::Row<size_t> label_vect_;        // sample labels
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
