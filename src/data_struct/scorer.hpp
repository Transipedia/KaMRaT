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
    kRSD3,
    kEntropy
};
const std::vector<std::string> kScorerNameVect{"ttest.padj", "ttest.pi", "SNR", "DIDS.score", "LR.acc", "Bayes.acc", "SVM.acc",
                                               "pearson", "spearman",
                                               "sd", "rsd1", "rsd2", "rsd3", "entropy"};

class Scorer
{
public:
    Scorer(const std::string &scorer_str, size_t nfold, const std::vector<std::string> &col_target_vect);

    const ScorerCode GetScorerCode() const;
    const std::string &GetScorerName() const;
    const double EstimateScore(std::vector<float> &count_vect) const;
    const double EstimateScore_old(const std::vector<float> &count_vect) const;

private:
    const ScorerCode scorer_code_;             // scoring method code
    const size_t nfold_;                       // prediction class number
    arma::Row<size_t> arma_categ_target_vect_; // armadillo categorical target vector
    std::vector<size_t> categ_target_vect_;
    std::vector<float> cntnu_target_vect_;     // continuous target vector
    size_t nclass_;                            // classification fold number

    const double LogTtestScore(const std::vector<float> &values) const;
    const double CalcSNRScore(const std::vector<float> & count_vect) const;
    const double CalcDIDSScore(const std::vector<float> count_vect) const;
    const double CalcPearsonScore(const std::vector<float> &count_vect) const;
    const double CalcSpearmanScore(const std::vector<float> &count_vect) const;
    const double CalcSDScore(const std::vector<float> & count_vect) const;
    const double CalcRSD1Score(const std::vector<float> &count_vect) const;
    const double CalcRSD2Score(const std::vector<float> &count_vect) const;
    /** WARNING: This method change the vector order to compute the median */
    const double CalcRSD3Score(std::vector<float> &count_vect) const;
    const double CalcEntropyScore(const std::vector<float> &count_vect) const;
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
