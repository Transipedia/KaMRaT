#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <ctime>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "data_struct/scorer.hpp"
#include "data_struct/count_tab.hpp"
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/rank.hpp"

void EvalScore(CountTab &feature_count_tab,
               seqVect_t &feature_vect,
               const std::string &count_tab_path,
               const std::string &smp_info_path,
               Scorer *scorer,
               const bool ln_transf,
               const bool standardize)
{
    std::ifstream count_tab_file(count_tab_path);
    if (!count_tab_file.is_open())
    {
        throw std::domain_error("count table file " + count_tab_path + " is not found");
    }
    size_t pos = count_tab_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && count_tab_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(count_tab_file);
    std::istream kmer_count_instream(&inbuf);

    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::string line;
    std::getline(kmer_count_instream, line);
    feature_count_tab.MakeSmpCond(smp_info_path);
    feature_count_tab.MakeColumnInfo(line, scorer->GetScoreCmd());
    std::vector<size_t> smp_labels;
    feature_count_tab.GetSmpLabels(smp_labels);
    scorer->LoadSampleLabel(smp_labels, feature_count_tab.GetNbCondition());
    // std::cout << "feature\tall.f1.binary\tall.f1.micro\tall.accuracy\tcv.f1.binary\tcv.f1.micro\tcv.accuracy\tcv.f1.na_to_0" << std::endl;
    //----- Dealing with Following k-mer Count Lines -----//
    const std::string idx_path = feature_count_tab.GetIndexPath();
    std::ofstream idx_file(idx_path); // create new index file because gz file does not support random access
    if (!idx_file.is_open())          // to ensure the file is opened
    {
        throw std::domain_error("error open file: " + idx_path);
    }
    for (uint64_t feature_serial(0); std::getline(kmer_count_instream, line); ++feature_serial)
    {
        std::istringstream conv(line);
        std::string feature_seq;
        conv >> feature_seq; // first column as feature (string)
        // std::cout << feature_seq << "\t";
        float feature_score;
        std::vector<float> count_vect;
        if (!feature_count_tab.AddRowAsString(feature_score, count_vect, line, idx_file) && scorer->GetScoreMethod() == "user") // firstly call AddRowAsString !
        {
            throw std::domain_error("user score column not found:" + scorer->GetScoreCmd());
        }
        else if (scorer->GetScoreMethod() != "user")
        {
            if (ln_transf)
            {
                for (size_t i(0); i < count_vect.size(); ++i)
                {
                    count_vect[i] = log(count_vect[i] + 1);
                }
            }
            feature_score = scorer->CalcScore(count_vect, standardize);
        }
        feature_vect.emplace_back(feature_seq, feature_serial, feature_score); // serial number as uniqcode
    }
    idx_file.close();
    count_tab_file.close();
}

void SortScore(seqVect_t &feature_vect, const std::string &sort_mode)
{
    if (sort_mode == "dec:abs")
    {
        auto comp = [](SeqElem feature1, SeqElem feature2) -> bool { return fabs(feature1.GetScore("origin")) > fabs(feature2.GetScore("origin")); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "dec")
    {
        auto comp = [](SeqElem feature1, SeqElem feature2) -> bool { return feature1.GetScore("origin") > feature2.GetScore("origin"); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc:abs")
    {
        auto comp = [](SeqElem feature1, SeqElem feature2) -> bool { return fabs(feature1.GetScore("origin")) < fabs(feature2.GetScore("origin")); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode == "inc")
    {
        auto comp = [](SeqElem feature1, SeqElem feature2) -> bool { return feature1.GetScore("origin") < feature2.GetScore("origin"); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else
    {
        throw std::domain_error("unknown sort mode: " + sort_mode);
    }
}

void PValueAdjustmentBH(seqVect_t &feature_vect)
{
    double tot_num = feature_vect.size();
    for (size_t i(0); i < tot_num; ++i)
    {
        const double adjust_factor = tot_num / (i + 1);
        feature_vect[i].ScaleScore(adjust_factor, 0, 1);
    }
    for (auto iter = feature_vect.rbegin() + 1; iter < feature_vect.rend(); ++iter)
    {
        double last_score = (iter - 1)->GetScore("final");
        if (iter->GetScore("final") > last_score)
        {
            iter->ScaleScore(1, last_score, last_score);
        }
    }
}

void ModelPrint(const seqVect_t &feature_vect,
                const std::string &idx_path,
                const size_t nb_sel,
                const CountTab &count_tab)
{
    std::cout << count_tab.GetColName(0) << "\tscore";
    for (size_t i(1); i < count_tab.GetNbColumn(); ++i)
    {
        std::cout << "\t" << count_tab.GetColName(i);
    }
    std::cout << std::endl;

    size_t parsed_nb_sel = (nb_sel == 0) ? feature_vect.size() : nb_sel;
    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("cannot open count index file: " + idx_path);
    }
    // std::cout << nb_sel << "\t" << parsed_nb_sel << std::endl;
    for (size_t i(0); i < parsed_nb_sel; ++i)
    {
        std::string row_string;
        count_tab.GetRowString(row_string, feature_vect[i].GetUniqCode(), idx_file, std::to_string(feature_vect[i].GetScore("final"))); // serial number was used as uniqcode
        std::cout << row_string << std::endl;
    }
    idx_file.close();
}

int RankMain(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string count_tab_path, smp_info_path, score_method, score_cmd, sort_mode, idx_path("./counts.idx");
    bool ln_transf(false), standardize(false);
    size_t nb_sel(0);

    ParseOptions(argc, argv, idx_path, smp_info_path, score_method, score_cmd, sort_mode, nb_sel, ln_transf, standardize, count_tab_path);
    Scorer *scorer;
    if (score_method.empty() || score_method == "sd")
    {
        scorer = new SDScorer(sort_mode);
    }
    else if (score_method == "rsd")
    {
        scorer = new RelatSDScorer(sort_mode);
    }
    else if (score_method == "ttest")
    {
        scorer = new TtestScorer(sort_mode);
    }
    else if (score_method == "es")
    {
        scorer = new EffectSizeScorer(sort_mode);
    }
    else if (score_method == "lfc")
    {
        if (score_cmd.empty())
        {
            scorer = new LFCScorer("mean", sort_mode);
        }
        else
        {
            scorer = new LFCScorer(score_cmd, sort_mode);
        }
    }
    else if (score_method == "nb")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = new NaiveBayesScorer(sort_mode, nb_fold);
    }
    else if (score_method == "rg")
    {
        size_t nb_fold = (score_cmd.empty() ? 1 : std::stoi(score_cmd));
        scorer = new RegressionScorer(sort_mode, nb_fold);
    }
    else if (score_method == "svm")
    {
        scorer = new SVMScorer(sort_mode);
    }
    else if (score_method == "user")
    {
        scorer = new UserScorer(sort_mode);
    }
    else
    {
        throw std::invalid_argument("unknown scoring method: " + score_method);
    }

    PrintRunInfo(count_tab_path, idx_path, smp_info_path, scorer->GetScoreMethod(), scorer->GetScoreCmd(), scorer->GetSortMode(), scorer->GetNbFold(),
                 nb_sel, ln_transf, standardize);
    CountTab count_tab(0, true, idx_path); // kamratRank does not care k-length and strandedness
    seqVect_t feature_vect;

    EvalScore(count_tab, feature_vect, count_tab_path, smp_info_path, scorer, ln_transf, standardize);
    SortScore(feature_vect, scorer->GetSortMode());
    if (scorer->GetScoreMethod() == "ttest")
    {
        PValueAdjustmentBH(feature_vect);
    }
    ModelPrint(feature_vect, idx_path, nb_sel, count_tab);

    delete scorer;

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
