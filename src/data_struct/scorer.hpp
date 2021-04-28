#ifndef KAMRAT_DATASTRUCT_SCORER_HPP
#define KAMRAT_DATASTRUCT_SCORER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <armadillo>

enum ScorerCode
{
    kTtest = 0,
    kSNR,
    kLR,
    kNBC,
    kSVM
};
const std::vector<std::string> kScorerNameVect{"padj.ttest", "SNR", "acc.LR", "acc.NBC", "acc.SVM"};

class Scorer
{
public:
    Scorer(const std::string &scorer_str, size_t nfold, const std::vector<double> &smp_sum_vect,
           const std::vector<size_t> &condi_vect, const std::vector<size_t> &batch_vect);

    const ScorerCode GetScorerCode() const;
    const std::string &GetScorerName() const;
    const double EstimateScore(const std::vector<float> &count_vect, bool no_norm, bool ln_transf, bool standardize) const;

private:
    const ScorerCode scorer_code_;      // scoring method code
    const size_t nfold_;                // prediction class number
    arma::Mat<double> arma_nf_vect_;    // normalization factor vector
    arma::Row<size_t> arma_condi_vect_; // real condition label vector
    size_t nclass_, nbatch_;            // classification fold number
    arma::Mat<double> arma_batch_vect_; // batch label vector
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
