#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <memory>
#include <ctime>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "data_struct/scorer.hpp"
#include "data_struct/tab_elem.hpp"
#include "data_struct/tab_header.hpp"
#include "data_struct/feature_elem.hpp"
#include "run_info_parser/rank.hpp"

void ScanCountComputeNF(featureVect_t &feature_vect, std::vector<double> &nf_vect, TabHeader &tab_header,
                        const std::string &raw_counts_path, const std::string &idx_path, const bool assign_score)
{
    std::ifstream raw_counts_file(raw_counts_path);
    if (!raw_counts_file.is_open())
    {
        throw std::domain_error("count table " + raw_counts_path + " was not found");
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = raw_counts_path.find_last_of(".");
        if (pos != std::string::npos && raw_counts_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(raw_counts_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line, str_x;
    std::getline(kmer_count_instream, line);
    tab_header.MakeColumnInfo(line, "");
    nf_vect.resize(tab_header.GetNbCount(), 0);
    std::cerr << "\t => Number of sample parsed: " << nf_vect.size() << std::endl;

    std::ofstream idx_file(idx_path);
    if (!idx_file.is_open()) // to ensure the file is opened
    {
        throw std::domain_error("error open file: " + idx_path);
    }

    std::vector<float> count_vect;
    std::istringstream conv;
    TabElem tab_elem;
    double mean_sample_sum(0); // mean of sample-sum vector
    while (std::getline(kmer_count_instream, line))
    {
        float score;
        conv.str(line);
        score = tab_elem.ParseTabElem(conv, idx_file, count_vect, tab_header);
        for (size_t i(0); i < count_vect.size(); ++i)
        {
            nf_vect[i] += (count_vect[i] + 1);      // + 1 for we will add an offset 1 to raw count for avoiding log(0)
            mean_sample_sum += (count_vect[i] + 1); // + 1 for we will add an offset 1 to raw count for avoiding log(0)
        }
        feature_vect.emplace_back(tab_elem.GetIndexPos(), score);
        count_vect.clear();
        conv.clear();
    }
    mean_sample_sum /= nf_vect.size(); // mean of sample-sum vector
    for (size_t i(0); i < nf_vect.size(); ++i)
    {
        nf_vect[i] /= mean_sample_sum;
    }
    raw_counts_file.close();
    idx_file.close();
}

void EvalScore(featureVect_t &feature_vect,
               std::ifstream &idx_file, const TabHeader &tab_header, const countTab_t &count_tab, const std::vector<double> nf_vect,
               const std::unique_ptr<Scorer> &scorer, const bool ln_transf, const bool standardize)
{
    std::vector<float> count_vect;
    for (size_t i_feature(0); i_feature < count_tab.size(); ++i_feature)
    {
        count_tab[i_feature].GetCountVect(count_vect, idx_file, tab_header.GetNbCount());
        for (size_t i_smp(0); i_smp < count_vect.size(); ++i_smp)
        {
            count_vect[i_smp] = (count_vect[i_smp] + 1) * nf_vect[i_smp];
        }
        scorer->LoadSampleCount(count_vect, ln_transf, standardize);
        if (scorer->GetScoreMethod() == "user") // user scoring mode but not find a score column
        {
            size_t rep_pos = tab_header.GetRepColPos();
            if (rep_pos == 0)
            {
                throw std::domain_error("user score column not found:" + scorer->GetScoreCmd());
            }
            feature_vect.emplace_back(i_feature,
                                      count_tab[i_feature].GetValueAt(idx_file, tab_header.GetColSerialAt(rep_pos), tab_header.GetNbCount()),
                                      scorer); // directly assign score with the given column
        }
        else
        {
            feature_vect.emplace_back(i_feature, scorer); // index position as uniqcode
        }
        count_vect.clear();
    }
}

void SortScore(featureVect_t &feature_vect, const std::string &sort_mode)
{
    if (sort_mode == "dec:abs")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return fabs(feature1.GetScore()) > fabs(feature2.GetScore()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "dec")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return feature1.GetScore() > feature2.GetScore(); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc:abs")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return fabs(feature1.GetScore()) < fabs(feature2.GetScore()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return feature1.GetScore() < feature2.GetScore(); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else
    {
        throw std::domain_error("unknown sort mode: " + sort_mode);
    }
}

void PValueAdjustmentBH(featureVect_t &feature_vect)
{
    double tot_num = feature_vect.size();
    for (size_t i(0); i < tot_num; ++i)
    {
        const double adjust_factor = tot_num / (i + 1);
        feature_vect[i].ScaleScore(adjust_factor, 0, 1);
    }
    for (auto iter = feature_vect.rbegin() + 1; iter < feature_vect.rend(); ++iter)
    {
        double last_score = (iter - 1)->GetScore();
        if (iter->GetScore() > last_score)
        {
            iter->ScaleScore(1, last_score, last_score);
        }
    }
}

void ModelPrint(const featureVect_t &feature_vect,
                const std::string &idx_path,
                const size_t nb_sel,
                const TabHeader &tab_header,
                const std::string &out_path)
{
    std::ofstream out_file;
    if (!out_path.empty())
    {
        out_file.open(out_path);
        if (!out_file.is_open())
        {
            throw std::domain_error("cannot open file: " + out_path);
        }
    }
    auto backup_buf = std::cout.rdbuf();
    if (!out_path.empty()) // output to file if a path is given, to screen if not
    {
        std::cout.rdbuf(out_file.rdbuf());
    }

    std::cout << tab_header.GetColNameAt(0) << "\tscore";
    for (size_t i(0); i < tab_header.GetNbCondition(); ++i)
    {
        std::cout << "\tnorm.mean" << static_cast<char>('A' + i);
    }
    for (size_t i(1); i < tab_header.GetNbCol(); ++i)
    {
        std::cout << "\t" << tab_header.GetColNameAt(i);
    }
    std::cout << std::endl;

    size_t parsed_nb_sel = (nb_sel == 0) ? feature_vect.size() : nb_sel;
    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("cannot open count index file: " + idx_path);
    }
    std::vector<float> count_vect, value_vect;
    for (size_t i(0); i < parsed_nb_sel; ++i)
    {
        size_t i_feature = feature_vect[i].GetSerialInTab();
        count_tab[i_feature].GetVectsAndClear(count_vect, value_vect, idx_file, tab_header.GetNbCount(), tab_header.GetNbValue());
        for (size_t j(1); j < tab_header.GetNbCol(); ++j)
        {
            size_t col_serial = tab_header.GetColSerialAt(j);
            if (tab_header.GetColNatureAt(i) >= 'A' && tab_header.GetColNatureAt(i) <= 'Z') // count column => output according to quant_mode
            {
                std::cout << "\t" << count_vect[col_serial];
            }
            else if (tab_header.GetColNatureAt(i) == 'v') // value column => output that related with rep-k-mer
            {
                std::cout << "\t" << value_vect[col_serial];
            }
        }
        std::cout << std::endl;
        std::cout << row_string.substr(0, row_string.find_first_of(" \t") + 1) << feature_vect[i].GetScore();
        for (const auto m : feature_vect[i].GetCondiMeans())
        {
            std::cout << "\t" << m;
        }
        std::cout << row_string.substr(row_string.find_first_of(" \t")) << std::endl;
    }

    idx_file.close();
    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
}

int RankMain(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string count_tab_path, smp_info_path, score_method, score_cmd, sort_mode, idx_path, nf_path, out_path;
    bool ln_transf(false), standardize(false);
    size_t nb_sel(0);

    ParseOptions(argc, argv, idx_path, nf_path, smp_info_path, score_method, score_cmd, sort_mode, nb_sel, ln_transf, standardize, out_path, count_tab_path);
    std::unique_ptr<Scorer> scorer;
    if (score_method.empty() || score_method == "sd")
    {
        scorer = std::make_unique<SDScorer>(sort_mode);
    }
    else if (score_method == "rsd")
    {
        scorer = std::make_unique<RelatSDScorer>(sort_mode);
    }
    else if (score_method == "ttest")
    {
        scorer = std::make_unique<TtestScorer>(sort_mode);
    }
    else if (score_method == "es")
    {
        scorer = std::make_unique<EffectSizeScorer>(sort_mode);
    }
    else if (score_method == "lfc")
    {
        if (score_cmd.empty())
        {
            scorer = std::make_unique<LFCScorer>("mean", sort_mode);
        }
        else
        {
            scorer = std::make_unique<LFCScorer>(score_cmd, sort_mode);
        }
    }
    else if (score_method == "nb")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = std::make_unique<NaiveBayesScorer>(sort_mode, nb_fold);
    }
    else if (score_method == "rg")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = std::make_unique<RegressionScorer>(sort_mode, nb_fold);
    }
    else if (score_method == "svm")
    {
        scorer = std::make_unique<SVMScorer>(sort_mode);
    }
    else if (score_method == "user")
    {
        scorer = std::make_unique<UserScorer>(sort_mode);
    }
    else
    {
        throw std::invalid_argument("unknown scoring method: " + score_method);
    }

    PrintRunInfo(count_tab_path, idx_path, nf_path, smp_info_path,
                 scorer->GetScoreMethod(), scorer->GetScoreCmd(), scorer->GetSortMode(), scorer->GetNbFold(),
                 nb_sel, ln_transf, standardize, out_path);

    TabHeader count_tab_header(smp_info_path);
    countTab_t count_tab;
    std::vector<double> sample_nf;
    IndexCountComputeNF(count_tab, sample_nf, count_tab_header, count_tab_path, idx_path);

    featureVect_t feature_vect;
    EvalScore(feature_vect, count_tab_header, count_tab, scorer, ln_transf, standardize, idx_path);
    SortScore(feature_vect, scorer->GetSortMode());
    if (score_method == "ttest")
    {
        PValueAdjustmentBH(feature_vect);
        std::cerr << "p-value adjusted by BH method..." << std::endl
                  << std::endl;
    }
    ModelPrint(feature_vect, idx_path, nb_sel, count_tab_header, out_path);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
