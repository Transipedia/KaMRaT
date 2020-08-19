#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <ctime>

#include "data_struct/scorer.hpp"
#include "data_struct/count_tab_by_string.hpp"
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/rank.hpp"

void EvalScore(CountTabByString &feature_count_tab,
               seqVect_t &feature_vect,
               const std::string &count_tab_path,
               const std::string &sample_info_path,
               Scorer *scorer)
{
    std::ifstream count_tab_file(count_tab_path);
    if (!count_tab_file.is_open())
    {
        throw std::domain_error("count table file " + count_tab_path + " is not found");
    }
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::string line;
    std::getline(count_tab_file, line);
    feature_count_tab.MakeSmpCond(sample_info_path);
    feature_count_tab.MakeColumnInfo(line, scorer->GetScoreCmd());
    std::vector<size_t> smp_labels;
    feature_count_tab.GetSmpLabels(smp_labels);
    scorer->LoadSampleLabel(smp_labels, feature_count_tab.GetNbCondition());
    //----- Dealing with Following k-mer Count Lines -----//
    for (size_t feature_serial(0), idx_pos(count_tab_file.tellg());
         std::getline(count_tab_file, line);
         ++feature_serial, idx_pos = count_tab_file.tellg())
    {
        float feature_score;
        std::vector<float> count_vect;
        if (!feature_count_tab.IndexWithString(feature_score, count_vect, line, idx_pos) && scorer->GetScoreMethod() == "user")
        {
            throw std::domain_error("user score column not found:" + scorer->GetScoreCmd());
        }
        else if (scorer->GetScoreMethod() != "user")
        {
            feature_score = scorer->CalcScore(count_vect);
        }
        std::istringstream conv(line);
        std::string feature_seq;
        conv >> feature_seq; // first column as feature (string)
        feature_vect.emplace_back(feature_seq, feature_serial, feature_score);
    }
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

void ModelPrint(const seqVect_t &feature_vect, const std::string &count_tab_path, const size_t top_num, const CountTabByString &count_tab)
{
    size_t parsed_top_num = (top_num == 0) ? feature_vect.size() : top_num;
    std::ifstream count_tab_file(count_tab_path);
    if (!count_tab_file.is_open())
    {
        throw std::domain_error("cannot open count index file: " + count_tab_path);
    }
    std::string line;
    std::getline(count_tab_file, line); // header line indicating sample names
    size_t first_delim = line.find_first_of("\t ");
    std::cout << line.substr(0, first_delim + 1) << "score" << line.substr(first_delim) << std::endl;
    for (size_t i(0); i < parsed_top_num; ++i)
    {
        count_tab.PrintFromInput(feature_vect[i].GetSerial(), count_tab_file, std::to_string(feature_vect[i].GetScore("final")));
    }
    count_tab_file.close();
}

int main(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string count_tab_path, sample_info_path, score_method, score_cmd, sort_mode;
    size_t nb_sel(0);

    ParseOptions(argc, argv, sample_info_path, score_method, score_cmd, sort_mode, nb_sel, count_tab_path);
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
        size_t nb_fold = (score_cmd.empty() ? 2 : std::stoi(score_cmd));
        scorer = new NaiveBayesScorer(sort_mode, nb_fold);
    }
    else if (score_method == "rg")
    {
        size_t nb_fold = (score_cmd.empty() ? 2 : std::stoi(score_cmd));
        scorer = new RegressionScorer(sort_mode, nb_fold);
    }
    else if (score_method == "user")
    {
        scorer = new UserScorer(sort_mode);
    }
    else
    {
        throw std::invalid_argument("unknown scoring method: " + score_method);
    }
    
    PrintRunInfo(count_tab_path, sample_info_path, scorer->GetScoreMethod(), scorer->GetScoreCmd(), scorer->GetSortMode(), scorer->GetNbFold(), nb_sel);
    CountTabByString count_tab;
    seqVect_t feature_vect;

    EvalScore(count_tab, feature_vect, count_tab_path, sample_info_path, scorer);
    SortScore(feature_vect, scorer->GetSortMode());
    if (scorer->GetScoreMethod() == "ttest")
    {
        PValueAdjustmentBH(feature_vect);
    }
    ModelPrint(feature_vect, count_tab_path, nb_sel, count_tab);

    delete scorer;

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
