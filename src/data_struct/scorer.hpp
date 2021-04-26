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
    Scorer(const std::string &scorer_str, size_t nfold, const std::vector<double> smp_sum_vect);

    const std::string &GetScorerName() const;
    const void LoadSampleLabels(const std::vector<size_t> &label_vect);
    const double EstimateScore(const std::vector<float> &count_vect, bool no_norm, bool ln_transf, bool standardize) const;

private:
    const ScorerCode scorer_code_;      // scoring method code
    const size_t nfold_;                // prediction class number
    arma::Row<size_t> arma_label_vect_; // real label vector
    size_t nclass_;                     // classification fold number
    arma::Mat<double> arma_nf_vect_;    // normalization factor vector
};

#endif //KAMRAT_DATASTRUCT_SCORER_HPP
