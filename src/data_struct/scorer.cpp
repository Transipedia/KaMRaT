#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include <mlpack/core/cv/k_fold_cv.hpp>
#include <mlpack/core/cv/metrics/accuracy.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
#include <mlpack/methods/linear_svm/linear_svm.hpp>

#include "scorer.hpp" // armadillo library need be included after mlpack

/* ============================== Scoring Method ========================= *\
 * ttest.padj    p-adj of t-test, adjusted by Benjamini-Hochberg procedure *
 * ttest.pi      pi-value of t-test                                        *
 * snr           signal-to-noise ratio                                     *
 * dids          DIDS score                                                *
 * lr            logistic regression                                       *
 * nbc           naive Bayes classifier                                    *
 * svm           support vector machine                                    *
 * sd            standard deviation (non-supervised)                       *
 * rsd           relative standard deviation (non-supervised)              *
\* ======================================================================= */

const ScorerCode ParseScorerCode(const std::string &scorer_str)
{
    if (scorer_str == "ttest.padj")
    {
        return ScorerCode::kTtestPadj;
    }
    if (scorer_str == "ttest.pi")
    {
        return ScorerCode::kTtestPi;
    }
    else if (scorer_str == "snr")
    {
        return ScorerCode::kSNR;
    }
    else if (scorer_str == "dids")
    {
        return ScorerCode::kDIDS;
    }
    else if (scorer_str == "lr")
    {
        return ScorerCode::kLR;
    }
    else if (scorer_str == "bayes")
    {
        return ScorerCode::kBayes;
    }
    else if (scorer_str == "svm")
    {
        return ScorerCode::kSVM;
    }
    else if (scorer_str == "sd")
    {
        return ScorerCode::kSD;
    }
    else if (scorer_str == "rsd")
    {
        return ScorerCode::kRSD;
    }
    else
    {
        throw std::invalid_argument("unknown ranking method: " + scorer_str);
    }
}

const double CalcTtestScore(const arma::Mat<double> &&arma_count_vect1, const arma::Mat<double> &&arma_count_vect2, const bool return_pi)
{
    double mean1 = arma::mean(arma::mean(arma_count_vect1)), mean2 = arma::mean(arma::mean(arma_count_vect2)),
           stddev1 = arma::mean(arma::stddev(arma_count_vect1)), stddev2 = arma::mean(arma::stddev(arma_count_vect2)),
           praw;
    size_t nb1 = arma_count_vect1.size(), nb2 = arma_count_vect2.size();
    if (stddev1 == 0 && stddev2 == 0)
    {
        praw = 1;
    }
    else
    {
        double t1 = stddev1 * stddev1 / nb1,
               t2 = stddev2 * stddev2 / nb2,
               df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (nb1 - 1) + t2 * t2 / (nb2 - 1)),
               t_stat = (mean1 - mean2) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        praw = 2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat)));
    }
    if (return_pi)
    {
        return (-log10(praw) * fabs(mean2 - mean1));
    }
    else
    {
        return praw;
    }
}

const double CalcSNRScore(const arma::Mat<double> &&arma_count_vect1, const arma::Mat<double> &&arma_count_vect2)
{
    double mean1 = arma::mean(arma::mean(arma_count_vect1)), mean2 = arma::mean(arma::mean(arma_count_vect2)),
           stddev1 = arma::mean(arma::stddev(arma_count_vect1)), stddev2 = arma::mean(arma::stddev(arma_count_vect2));
    if (stddev1 == 0 && stddev2 == 0)
    {
        return 0;
    }
    else
    {
        return ((mean1 - mean2) / (stddev1 + stddev2));
    }
}

const double CalcDIDSScore(const arma::Row<size_t> &arma_label_vect, const arma::Mat<double> &arma_count_vect, const size_t nclass)
{
    double max_score(0), sqrt_sum, ref_max;
    for (size_t i_condi(0); i_condi < nclass; ++i_condi)
    {
        ref_max = arma_count_vect.elem(arma::find(arma_label_vect == i_condi)).max();
        sqrt_sum = 0;
        for (double x : arma::conv_to<std::vector<double>>::from(arma_count_vect.elem(arma::find(arma_label_vect != i_condi))))
        {
            sqrt_sum += (x > ref_max ? sqrt(x - ref_max) : 0); // alternatively, sqrt((x - ref_max) / ref_max)
        }
        max_score = (sqrt_sum > max_score ? sqrt_sum : max_score);
    }
    return max_score;
}

const double CalcLRScore(const size_t nfold, const arma::Row<size_t> &arma_label_vect, const arma::Mat<double> &arma_count_vect)
{
    if (nfold == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::regression::LogisticRegression<> lr(arma_count_vect, arma_label_vect);
        mlpack::cv::Accuracy acc;
        return acc.Evaluate(lr, arma_count_vect, arma_label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>, mlpack::cv::Accuracy>
            score_data(nfold, arma_count_vect, arma_label_vect);
        return score_data.Evaluate();
    }
}

const double CalcBayesScore(const size_t nfold, const arma::Row<size_t> &arma_label_vect, const arma::Mat<double> &arma_count_vect, const size_t nclass)
{
    if (nfold == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(arma_count_vect, arma_label_vect, nclass);
        mlpack::cv::Accuracy acc;
        return acc.Evaluate(nbc, arma_count_vect, arma_label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, mlpack::cv::Accuracy>
            score_data(nfold, arma_count_vect, arma_label_vect, nclass);
        return score_data.Evaluate();
    }
}

const double CalcSVMScore(const size_t nfold, const arma::Row<size_t> &arma_label_vect, const arma::Mat<double> &arma_count_vect, const size_t nclass)
{
    if (nfold == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::svm::LinearSVM<> lsvm(arma_count_vect, arma_label_vect, nclass);
        mlpack::cv::Accuracy acc;
        return acc.Evaluate(lsvm, arma_count_vect, arma_label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::svm::LinearSVM<>, mlpack::cv::Accuracy>
            score_data(nfold, arma_count_vect, arma_label_vect, nclass);
        return score_data.Evaluate();
    }
}

const double CalcSDScore(const arma::Mat<double> &arma_count_vect)
{
    return arma::mean(arma::stddev(arma_count_vect, 0, 1)); // arma_count_vect is not processed by .elem(), so still row vectors
}

const double CalcRSDScore(const arma::Mat<double> &arma_count_vect)
{
    double sd = arma::mean(arma::stddev(arma_count_vect, 0, 1)),
           mean = arma::mean(arma::mean(arma_count_vect, 1)); // arma_count_vect is not processed by .elem(), so still row vectors
    return (mean == 0 ? 0 : (sd / mean));
}

Scorer::Scorer(const std::string &scorer_str, const size_t nfold,
               const std::vector<size_t> &condi_label_vect, const std::vector<size_t> &batch_label_vect)
    : scorer_code_(ParseScorerCode(scorer_str)), nfold_((nfold == 0 ? condi_label_vect.size() : nfold)), nbatch_(0)
{
    if (!condi_label_vect.empty())
    {
        arma_condi_vect_ = arma::conv_to<arma::Row<size_t>>::from(condi_label_vect);
        // arma_condi_vect_.print("Label vector:");
        nclass_ = arma_condi_vect_.max() + 1;
    }
    if (nclass_ != 2 && (scorer_code_ == ScorerCode::kTtestPadj || scorer_code_ == ScorerCode::kTtestPi ||
                         scorer_code_ == ScorerCode::kSNR || scorer_code_ == ScorerCode::kLR))
    {
        throw std::domain_error("scoring by t-test, SNR, and LR only accept binary sample condition: " + std::to_string(nclass_));
    }
    if (nclass_ < 2 && (scorer_code_ == ScorerCode::kDIDS || scorer_code_ == ScorerCode::kBayes || scorer_code_ == ScorerCode::kSVM))
    {
        throw std::domain_error("scoring by DIDS, Bayes or SVM only accepts condition number >= 2");
    }
    if (!batch_label_vect.empty())
    {
        arma_batch_vect_ = arma::conv_to<arma::Row<double>>::from(batch_label_vect);
        // arma_batch_vect_.print("Batch vector:");
        nbatch_ = arma_batch_vect_.max() + 1;
    }
    if (nbatch_ > 1 && (scorer_code_ == ScorerCode::kTtestPadj || scorer_code_ == ScorerCode::kTtestPi ||
                        scorer_code_ == ScorerCode::kSNR || scorer_code_ == ScorerCode::kDIDS || scorer_code_ == ScorerCode::kSD))
    {
        throw std::invalid_argument("T-test, SNR and DIDS do not support batch effect correction");
    }
}

const ScorerCode Scorer::GetScorerCode() const
{
    return scorer_code_;
}

const std::string &Scorer::GetScorerName() const
{
    return kScorerNameVect[scorer_code_];
}

const double Scorer::EstimateScore(const std::vector<float> &count_vect) const
{
    static arma::Mat<double> arma_count_vect;
    arma_count_vect = arma::conv_to<arma::Row<double>>::from(count_vect);
    // arma_count_vect.print("Count vector before transformation: ");
    if (scorer_code_ == ScorerCode::kTtestPadj || scorer_code_ == ScorerCode::kTtestPi) // if t-test, apply log2(x + 1) transformation
    {
        arma_count_vect = log2(arma_count_vect + 1);
    }
    else if (scorer_code_ == ScorerCode::kLR || scorer_code_ == ScorerCode::kBayes || scorer_code_ == ScorerCode::kSVM) // if ML, apply standarization
    {
        arma_count_vect = (arma_count_vect - arma::mean(arma::mean(arma_count_vect, 1))) / arma::mean(arma::stddev(arma_count_vect, 0, 1));
    }
    if (scorer_code_ == ScorerCode::kLR || scorer_code_ == ScorerCode::kBayes || scorer_code_ == ScorerCode::kSVM) // add batch rows for LR, Bayes, SVM
    {
        arma_count_vect = arma::join_cols(arma_count_vect, arma_batch_vect_);
    }
    // arma_count_vect.print("Count vector after transformation:");

    switch (scorer_code_)
    {
    case ScorerCode::kTtestPadj:
    case ScorerCode::kTtestPi:
        return CalcTtestScore(arma_count_vect.elem(arma::find(arma_condi_vect_ == 0)),
                              arma_count_vect.elem(arma::find(arma_condi_vect_ == 1)),
                              scorer_code_ == ScorerCode::kTtestPi);
    case ScorerCode::kSNR:
        return CalcSNRScore(arma_count_vect.elem(arma::find(arma_condi_vect_ == 0)),
                            arma_count_vect.elem(arma::find(arma_condi_vect_ == 1)));
    case ScorerCode::kDIDS:
        return CalcDIDSScore(arma_condi_vect_, arma_count_vect, nclass_);
    case ScorerCode::kLR:
        return CalcLRScore(nfold_, arma_condi_vect_, arma_count_vect);
    case ScorerCode::kBayes:
        return CalcBayesScore(nfold_, arma_condi_vect_, arma_count_vect, nclass_);
    case ScorerCode::kSVM:
        return CalcSVMScore(nfold_, arma_condi_vect_, arma_count_vect, nclass_);
    case ScorerCode::kSD:
        return CalcSDScore(arma_count_vect);
    case ScorerCode::kRSD:
        return CalcRSDScore(arma_count_vect);
    default:
        return std::nan("");
    }
}
