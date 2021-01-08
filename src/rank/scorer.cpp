#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include <mlpack/core/cv/k_fold_cv.hpp>
#include <mlpack/core/cv/metrics/f1.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/methods/linear_svm/linear_svm.hpp>
#include <mlpack/methods/linear_svm/linear_svm_function.hpp>

#include "scorer.hpp" // armadillo library need be included after mlpack

/* ============================ Scoring Method ============================ *\
 * rsd          relative standard deviation                                  *
 * ttest        p-value of t-test, adjusted by Benjamini-Hochberg procedure  *
 * snr          signal-to-noise ratio                                        *
 * lrc          logistic regression classifier                               *
 * nbc          naive Bayes classifier                                       *
 * svm          support vector machine                                       *
 * user:name    use representative value in table as score                   *
\* ======================================================================== */

const SortModeCode ParseSortCodeFromScoreCode(const ScoreMethodCode score_code)
{
    switch (score_code)
    {
    case ScoreMethodCode::kRelatSD:
    case ScoreMethodCode::kNaiveBayes:
    case ScoreMethodCode::kLogitReg:
        return SortModeCode::kDec;
    case ScoreMethodCode::kSNR:
        return SortModeCode::kDecAbs;
    case ScoreMethodCode::kTtest:
    case ScoreMethodCode::kSVM:
        return SortModeCode::kInc;
    default:
        return SortModeCode::kDec;
    }
}

const double CalcRelatSDScore(const double all_mean, const double all_stddev)
{
    if (all_mean <= 1)
    {
        return all_stddev;
    }
    else
    {
        return (all_stddev / all_mean);
    }
}

const double CalcTtestScore(const size_t nb1, const size_t nb2,
                            const double mean1, const double mean2, const double stddev1, const double stddev2)
{
    if (stddev1 == 0 && stddev2 == 0)
    {
        return 1;
    }
    else
    {
        double t1 = stddev1 * stddev1 / nb1,
               t2 = stddev2 * stddev2 / nb2,
               df = (t1 + t2) * (t1 + t2) / (t1 * t1 / (nb1 - 1) + t2 * t2 / (nb2 - 1)),
               t_stat = (mean1 - mean2) / sqrt(t1 + t2);
        boost::math::students_t dist(df);
        return (2 * boost::math::cdf(boost::math::complement(dist, fabs(t_stat))));
    }
}

const double CalcSNRScore(const double mean1, const double mean2, const double stddev1, const double stddev2)
{
    if (stddev1 == 0 && stddev2 == 0)
    {
        return std::nan("");
    }
    else
    {
        return ((mean1 - mean2) / (stddev1 + stddev2));
    }
}

class MLMetrics
{
public:
    template <typename MLAlgorithm, typename DataType>
    static double Evaluate(MLAlgorithm &model, const DataType &data, const arma::Row<size_t> &labels)
    {
        if (labels.max() == 1) // binary conditions
        {
            // data.print("Data:");
            double score = mlpack::cv::F1<mlpack::cv::Binary>().Evaluate(model, data, labels);
            return (isnan(score) ? 0 : score);
        }
        else // multiple conditions
        {
            double score = mlpack::cv::F1<mlpack::cv::Micro>().Evaluate(model, data, labels);
            return score;
        }
    }
};

const double CalcLRCScore(const size_t nb_fold, const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect)
{
    const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::regression::LogisticRegression<> rgc(norm_count_vect, label_vect);
        MLMetrics my_f1;
        return my_f1.Evaluate(rgc, norm_count_vect, label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::regression::LogisticRegression<>, MLMetrics>
            score_data(nb_fold_final, norm_count_vect, label_vect);
        return score_data.Evaluate();
    }
}

const double CalcNBCScore(const size_t nb_fold, const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect, const size_t nb_class)
{
    const size_t nb_fold_final = (nb_fold == 0 ? label_vect.size() : nb_fold);
    if (nb_fold_final == 1) // without cross-validation, train and test on the whole set
    {
        mlpack::naive_bayes::NaiveBayesClassifier<> nbc(norm_count_vect, label_vect, nb_class);
        MLMetrics my_f1;
        return my_f1.Evaluate(nbc, norm_count_vect, label_vect);
    }
    else // k-fold cross-validation (k=0 for leave-one-out cross-validation)
    {
        mlpack::cv::KFoldCV<mlpack::naive_bayes::NaiveBayesClassifier<>, MLMetrics>
            score_data(nb_fold_final, norm_count_vect, label_vect, nb_class);
        return score_data.Evaluate();
    }
}

const double CalcSVMScore(const arma::Row<size_t> &label_vect, const arma::Mat<double> &norm_count_vect, const size_t nb_class)
{
    mlpack::svm::LinearSVM<> lsvm(norm_count_vect, label_vect, nb_class);
    mlpack::svm::LinearSVMFunction<> lsvm_fun(norm_count_vect, label_vect, nb_class);
    return lsvm_fun.Evaluate(lsvm.Parameters());
}

Scorer::Scorer(const ScoreMethodCode score_method_code) // for relat.sd, t-test, snr, svm
    : score_method_code_(score_method_code), rep_colname_(""), sort_mode_code_(ParseSortCodeFromScoreCode(score_method_code)), nb_fold_(0)
{
}

Scorer::Scorer(const ScoreMethodCode score_method_code, const size_t nb_fold) // for naive bayes, logistic regression
    : score_method_code_(score_method_code), rep_colname_(""), sort_mode_code_(ParseSortCodeFromScoreCode(score_method_code)), nb_fold_(nb_fold)
{
}

Scorer::Scorer(const std::string &rep_colname, const SortModeCode sort_mode_code) // column name given by user
    : score_method_code_(ScoreMethodCode::kUser), rep_colname_(rep_colname), sort_mode_code_(sort_mode_code), nb_fold_(0)
{
}

const void Scorer::LoadSampleLabels(const std::vector<size_t> &label_vect)
{
    label_vect_ = arma::conv_to<arma::Row<size_t>>::from(label_vect);
    // label_vect_.print("Label vector:");
    nb_class_ = label_vect_.max() + 1;
    if (score_method_code_ == ScoreMethodCode::kTtest) // only t-test needs
    {
        nbsmp_condi_.resize(nb_class_, 0);
        for (size_t x : label_vect)
        {
            nbsmp_condi_[x]++;
        }
    }
    if (nb_class_ != 2 &&
        (score_method_code_ == ScoreMethodCode::kTtest || score_method_code_ == ScoreMethodCode::kSNR || score_method_code_ == ScoreMethodCode::kLogitReg))
    {
        throw std::domain_error("scoring by T-test, signal-to-noise ratio, and logistic regression only accept binary sample condition");
    }
    if (nb_class_ < 2 &&
        (score_method_code_ == ScoreMethodCode::kNaiveBayes || score_method_code_ == ScoreMethodCode::kSVM))
    {
        throw std::domain_error("scoring by naive Bayes or SVM only accepts condition number >= 2");
    }
}

const ScoreMethodCode Scorer::GetScoreMethodCode() const
{
    return score_method_code_;
}

const std::string &Scorer::GetRepColname() const
{
    return rep_colname_;
}

const SortModeCode Scorer::GetSortModeCode() const
{
    return sort_mode_code_;
}

const size_t Scorer::GetNbFold() const
{
    return nb_fold_;
}

const size_t Scorer::GetNbClass() const
{
    return nb_class_;
}

const void Scorer::PrepareCountVect(const FeatureElem &feature_elem, const std::vector<double> &nf_vect, std::ifstream &idx_file,
                                    const bool to_ln, const bool to_standardize, const bool no_norm)
{
    static std::vector<float> raw_count_vect;
    const size_t kNbCount = label_vect_.size();
    norm_count_vect_.zeros(1, kNbCount);
    feature_elem.RetrieveCountVect(raw_count_vect, idx_file, kNbCount);
    if (no_norm) // first normalization
    {
        for (size_t i(0); i < kNbCount; ++i)
        {
            norm_count_vect_(0, i) = raw_count_vect[i];
        }
    }
    else
    {
        for (size_t i(0); i < kNbCount; ++i)
        {
            norm_count_vect_(0, i) = nf_vect[i] * raw_count_vect[i];
        }
    }
    if (to_ln) // then ln(x + 1) transformation
    {
        for (size_t i(0); i < kNbCount; ++i)
        {
            norm_count_vect_(0, i) = log(nf_vect[i] * raw_count_vect[i] + 1);
        }
    }
    if (to_standardize) // finally standardization
    {
        const double all_mean = arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect_)),
                     all_stddev = arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect_));
        for (size_t i(0); i < kNbCount; ++i)
        {
            norm_count_vect_(0, i) = (norm_count_vect_(0, i) - all_mean) / all_stddev;
        }
    }
    // norm_count_vect_.print("Norm count vector:");
}

const void Scorer::CalcFeatureStats(FeatureElem &feature_elem)
{
    if (score_method_code_ == ScoreMethodCode::kRelatSD)
    {
        feature_elem.ReserveCondiStats(1);
        const double all_mean = arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect_)),
                     all_stddev = arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect_));
        feature_elem.AddCondiStats(all_mean, all_stddev);
    }
    else
    {
        feature_elem.ReserveCondiStats(nb_class_);
        for (size_t i_class(0); i_class < nb_class_; ++i_class)
        {
            arma::Mat<double> &&norm_count_vect_i = norm_count_vect_.elem(arma::find(label_vect_ == i_class));
            const double condi_mean = arma::mean(arma::conv_to<arma::vec>::from(norm_count_vect_i)),
                         condi_stddev = arma::stddev(arma::conv_to<arma::vec>::from(norm_count_vect_i));
            feature_elem.AddCondiStats(condi_mean, condi_stddev);
        }
    }
}

const double Scorer::EvaluateScore(const FeatureElem &feature_elem) const
{
    switch (score_method_code_)
    {
    case ScoreMethodCode::kRelatSD:
        return CalcRelatSDScore(feature_elem.GetCondiMeanAt(0), feature_elem.GetCondiStddevAt(0));
    case ScoreMethodCode::kTtest:
        return CalcTtestScore(nbsmp_condi_[0], nbsmp_condi_[1],
                              feature_elem.GetCondiMeanAt(0), feature_elem.GetCondiMeanAt(1),
                              feature_elem.GetCondiStddevAt(0), feature_elem.GetCondiStddevAt(1));
    case ScoreMethodCode::kSNR:
        return CalcSNRScore(feature_elem.GetCondiMeanAt(0), feature_elem.GetCondiMeanAt(1),
                            feature_elem.GetCondiStddevAt(0), feature_elem.GetCondiStddevAt(1));
    case ScoreMethodCode::kLogitReg:
        return CalcLRCScore(nb_fold_, label_vect_, norm_count_vect_);
    case ScoreMethodCode::kNaiveBayes:
        return CalcNBCScore(nb_fold_, label_vect_, norm_count_vect_, nb_class_);
    case ScoreMethodCode::kSVM:
        return CalcSVMScore(label_vect_, norm_count_vect_, nb_class_);
    default:
        return std::nan("");
    }
}
