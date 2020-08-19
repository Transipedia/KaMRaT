#ifndef KAMRAT_UTILS_SCORER_HPP
#define KAMRAT_UTILS_SCORER_HPP

#include <string>
#include <vector>
#include <cmath>

#include "mlpack/core/cv/k_fold_cv.hpp"
#include "mlpack/core/cv/metrics/f1.hpp"
#include "mlpack/methods/naive_bayes/naive_bayes_classifier.hpp"
#include "mlpack/methods/softmax_regression/softmax_regression.hpp"
#include "armadillo"
#include "boost/math/distributions/students_t.hpp"

class Scorer
{
public:
    Scorer(const std::string &score_method, const std::string &score_cmd, const std::string &sort_mode);
    void LoadSampleLabel(const std::vector<size_t> &label_vect, const size_t nb_class);
    const std::string &GetScoreMethod() const;
    const std::string &GetSortMode() const;
    const std::string &GetScoreCmd() const;
    const size_t GetNbFold() const;
    const size_t GetNbClass() const;
    const float CalcScore(const std::vector<float> &sample_counts) const;

private:
    const std::string score_method_;  // nb, lr, sd, rsd, ttest, es, lfc, user
    const std::string sort_mode_;     // dec, dec::abs, inc, inc::abs
    const std::string score_cmd_;     // to nb_fold_ if nb or lr; else for mean or median in lfc or for score colname in user
    const size_t nb_fold_;            // for naive Bayes and logistic regression
    size_t nb_class_;                 // number of conditions
    arma::Row<size_t> sample_labels_; // sample labels
};

#endif //KAMRAT_UTILS_SCORER_HPP
