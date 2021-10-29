#ifndef KAMRAT_DATASTRUCT_SCORER_HPP
#define KAMRAT_DATASTRUCT_SCORER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <armadillo>

enum ScorerCode
{
    kTtestPadj = 0,
    kTtestPi,
    kSNR,
    kDIDS,
    kLR,
    kBayes,
    kSVM,
    kPearson,
    kSpearman,
    kSD,
    kRSD,
    kEntropy
};
const std::vector<std::string> kScorerNameVect{"ttest.padj", "ttest.pi", "SNR", "DIDS.score", "LR.acc", "Bayes.acc", "SVM.acc",
                                               "pearson", "spearman",
                                               "stddev", "relat.stddev", "entropy"};

class Scorer
{
public:
    Scorer(const std::string &scorer_str, size_t nfold, const std::vector<std::string> &col_target_vect, const std::vector<std::string> &col_batch_vect);
    // Scorer(const std::string &scorer_str, size_t nfold, const std::vector<size_t> &condi_vect, const std::vector<size_t> &batch_vect);

    const ScorerCode GetScorerCode() const;
    const std::string &GetScorerName() const;
    const double EstimateScore(const std::vector<float> &count_vect) const;

private:
    const ScorerCode scorer_code_;             // scoring method code
    const size_t nfold_;                       // prediction class number
    arma::Row<size_t> arma_categ_target_vect_; // armadillo categorical target vector
    std::vector<float> cntnu_target_vect_;     // continuous target vector
    size_t nclass_, nbatch_;                   // classification fold number
    arma::Mat<double> arma_batch_vect_;        // batch label vector
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
