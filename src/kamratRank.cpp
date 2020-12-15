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
#include "data_struct/tab_header.hpp"
#include "data_struct/feature_elem.hpp"
#include "run_info_parser/rank.hpp"

const void ScanCountComputeNF(featureVect_t &feature_vect, std::vector<double> &nf_vect, TabHeader &tab_header,
                              const std::string &raw_counts_path, const std::string &idx_path)
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
    std::istringstream conv(line);
    tab_header.MakeColumnInfo(conv, "");
    conv.clear();

    std::ofstream idx_file(idx_path);
    if (!idx_file.is_open()) // to ensure the file is opened
    {
        throw std::domain_error("error open file: " + idx_path);
    }
    
    std::vector<float> count_vect;
    nf_vect.resize(tab_header.GetNbCount(), 0);
    std::cerr << "\t => Number of sample parsed: " << nf_vect.size() << std::endl;

    double mean_sample_sum(0); // mean of sample-sum vector
    while (std::getline(kmer_count_instream, line))
    {
        conv.str(line);
        feature_vect.emplace_back(conv, count_vect, idx_file, tab_header);
        for (size_t i(0); i < count_vect.size(); ++i)
        {
            nf_vect[i] += (count_vect[i] + 1);      // + 1 for we will add an offset 1 to raw count for avoiding log(0)
            mean_sample_sum += (count_vect[i] + 1); // + 1 for we will add an offset 1 to raw count for avoiding log(0)
        }
        count_vect.clear();
        conv.clear();
    }
    mean_sample_sum /= nf_vect.size(); // mean of sample-sum vector
    for (size_t i(0); i < nf_vect.size(); ++i)
    {
        nf_vect[i] = mean_sample_sum / nf_vect[i];
    }
    raw_counts_file.close();
    idx_file.close();
}

void PrintNF(const std::string &smp_sum_outpath,
             const std::vector<double> &nf_vect,
             const TabHeader &tab_header)
{
    std::ofstream sum_out(smp_sum_outpath);
    for (size_t i(1), j(0); i < tab_header.GetNbCol(); ++i)
    {
        if (tab_header.IsCount(i))
        {
            sum_out << tab_header.GetColNameAt(i) << "\t" << nf_vect[j++] << std::endl;
        }
    }
    sum_out.close();
}

void EvalScore(featureVect_t &feature_vect,
               std::ifstream &idx_file, const std::vector<double> nf_vect,
               std::unique_ptr<Scorer> &scorer, const TabHeader tab_header,
               const bool ln_transf, const bool standardize)
{
    size_t nb_count = nf_vect.size();
    std::vector<float> count_vect;
    scorer->LoadSampleLabel(tab_header);
    for (size_t i_feature(0); i_feature < feature_vect.size(); ++i_feature)
    {
        feature_vect[i_feature].GetCountVect(count_vect, idx_file, nb_count);
        scorer->LoadSampleCount(count_vect, nf_vect);
        if (scorer->GetScoreMethod() != "user") // user scoring mode but not find a score column
        {
            feature_vect[i_feature].EvalFeatureElem(scorer);
        }
        count_vect.clear();
    }
}

void SortScore(featureVect_t &feature_vect, const std::string &sort_mode)
{
    if (sort_mode == "dec:abs")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return fabs(feature1.GetValue()) > fabs(feature2.GetValue()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "dec")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return feature1.GetValue() > feature2.GetValue(); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc:abs")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return fabs(feature1.GetValue()) < fabs(feature2.GetValue()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc")
    {
        auto comp = [](FeatureElem feature1, FeatureElem feature2) -> bool { return feature1.GetValue() < feature2.GetValue(); };
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
        feature_vect[i].ScaleValue(adjust_factor, 0, 1);
    }
    for (auto iter = feature_vect.rbegin() + 1; iter < feature_vect.rend(); ++iter)
    {
        double last_score = (iter - 1)->GetValue();
        if (iter->GetValue() > last_score)
        {
            iter->ScaleValue(1, last_score, last_score);
        }
    }
}

void ModelPrint(featureVect_t &feature_vect,
                std::ifstream &idx_file,
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
    std::string row_str;
    tab_header.MakeOutputHeaderStr(row_str);
    std::cout << row_str.substr(row_str.find_first_of("\t")) << std::endl;

    size_t parsed_nb_sel = (nb_sel == 0) ? feature_vect.size() : nb_sel;
    std::string out_row_str;
    size_t split_pos;
    for (size_t i(0); i < parsed_nb_sel; ++i)
    {
        feature_vect[i].MakeOutputRowStr(out_row_str, idx_file, tab_header.GetNbCount());
        split_pos = out_row_str.find_first_of("\t");
        std::cout << out_row_str.substr(0, split_pos + 1) << feature_vect[i].GetValue();
        for (const float m : feature_vect[i].GetNormCondiMeans())
        {
            std::cout << "\t" << m;
        }
        std::cout << out_row_str.substr(split_pos) << std::endl;
    }

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
}

int RankMain(int argc, char *argv[])
{
    std::clock_t begin_time = clock(), inter_time;
    std::string count_tab_path, smp_info_path, score_method, score_cmd, sort_mode, idx_path, nf_path, out_path;
    bool ln_transf(false), standardize(false);
    size_t nb_sel(0);

    ParseOptions(argc, argv, idx_path, nf_path, smp_info_path, score_method, score_cmd, sort_mode, nb_sel, ln_transf, standardize, out_path, count_tab_path);
    std::unique_ptr<Scorer> scorer;
    if (score_method.empty() || score_method == "sd")
    {
        scorer = std::make_unique<SDScorer>(sort_mode, ln_transf, standardize);
    }
    else if (score_method == "rsd")
    {
        scorer = std::make_unique<RelatSDScorer>(sort_mode, ln_transf, standardize);
    }
    else if (score_method == "ttest")
    {
        scorer = std::make_unique<TtestScorer>(sort_mode, ln_transf, standardize);
    }
    else if (score_method == "es")
    {
        scorer = std::make_unique<EffectSizeScorer>(sort_mode, ln_transf, standardize);
    }
    else if (score_method == "lfc")
    {
        if (score_cmd.empty())
        {
            scorer = std::make_unique<LFCScorer>("mean", sort_mode, ln_transf, standardize);
        }
        else
        {
            scorer = std::make_unique<LFCScorer>(score_cmd, sort_mode, ln_transf, standardize);
        }
    }
    else if (score_method == "nb")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = std::make_unique<NaiveBayesScorer>(sort_mode, nb_fold, ln_transf, standardize);
    }
    else if (score_method == "rg")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = std::make_unique<RegressionScorer>(sort_mode, nb_fold, ln_transf, standardize);
    }
    else if (score_method == "svm")
    {
        scorer = std::make_unique<SVMScorer>(sort_mode, ln_transf, standardize);
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

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    TabHeader count_tab_header(smp_info_path);
    featureVect_t feature_vect;
    std::vector<double> smp_nf_vect;
    ScanCountComputeNF(feature_vect, smp_nf_vect, count_tab_header, count_tab_path, idx_path);
    PrintNF(nf_path, smp_nf_vect, count_tab_header);
    std::cerr << "Count table scanning finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("error open sample count index file: " + idx_path);
    }
    EvalScore(feature_vect, idx_file, smp_nf_vect, scorer, count_tab_header, ln_transf, standardize);
    std::cerr << "Score evalution finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    SortScore(feature_vect, scorer->GetSortMode());
    std::cerr << "Score sorting finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    if (score_method == "ttest")
    {
        PValueAdjustmentBH(feature_vect);
        std::cerr << "p-value adjusted by BH method..." << std::endl
                  << std::endl;
        std::cerr << "P-value adjusting finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
        inter_time = clock();
    }

    ModelPrint(feature_vect, idx_file, nb_sel, count_tab_header, out_path);
    std::cerr << "Output finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    idx_file.close();
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
