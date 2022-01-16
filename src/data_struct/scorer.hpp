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
    kRSD1,
    kRSD2,
    kEntropy
};
const std::vector<std::string> kScorerNameVect{"ttest.padj", "ttest.pi", "SNR", "DIDS.score", "LR.acc", "Bayes.acc", "SVM.acc",
                                               "pearson", "spearman",
                                               "sd", "rsd1", "rsd2", "entropy"};

class Scorer
{
public:
    Scorer(const std::string &scorer_str, size_t nfold, const std::vector<std::string> &col_target_vect);

    const ScorerCode GetScorerCode() const;
    const std::string &GetScorerName() const;
    const double EstimateScore(const std::vector<float> &count_vect) const;

private:
    const ScorerCode scorer_code_;             // scoring method code
    const size_t nfold_;                       // prediction class number
    arma::Row<size_t> arma_categ_target_vect_; // armadillo categorical target vector
    std::vector<float> cntnu_target_vect_;     // continuous target vector
    size_t nclass_;                            // classification fold number
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
