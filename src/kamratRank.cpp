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
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/rank.hpp"

void EvalScore(seqVect_t &feature_vect,
               TabHeader &tab_header,
               const std::string &count_tab_path,
               const std::unique_ptr<Scorer> &scorer,
               const bool ln_transf,
               const bool standardize,
               const std::string &idx_path)
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
    tab_header.MakeColumnInfo(line, (scorer->GetScoreMethod() == "user" ? scorer->GetScoreCmd() : ""));
    scorer->LoadSampleLabel(tab_header);
    //----- Dealing with Following k-mer Count Lines -----//
    if (idx_path.empty())
    {
        throw std::domain_error("index file path not provided");
    }
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
        TabElem tab_elem(conv, idx_file, count_vect, feature_score, tab_header);
        if (tab_header.GetRepColPos() == 0 && scorer->GetScoreMethod() == "user") // user scoring mode but not find a score column
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
        feature_vect.emplace_back(feature_seq, tab_elem.GetIndexPos(), feature_score); // index position as uniqcode
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
    // std::cout << nb_sel << "\t" << parsed_nb_sel << std::endl;
    for (size_t i(0); i < parsed_nb_sel; ++i)
    {
        std::string row_string;
        idx_file.seekg(feature_vect[i].GetUniqCode());
        std::getline(idx_file, row_string);
        std::cout << row_string.substr(0, row_string.find_first_of(" \t") + 1)
                  << feature_vect[i].GetScore("final")
                  << row_string.substr(row_string.find_first_of(" \t")) << std::endl;
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
    std::string count_tab_path, smp_info_path, score_method, score_cmd, sort_mode, idx_path("./counts.idx"), out_path;
    bool ln_transf(false), standardize(false);
    size_t nb_sel(0);

    ParseOptions(argc, argv, idx_path, smp_info_path, score_method, score_cmd, sort_mode, nb_sel, ln_transf, standardize, out_path, count_tab_path);
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

    PrintRunInfo(count_tab_path, idx_path, smp_info_path, scorer->GetScoreMethod(), scorer->GetScoreCmd(), scorer->GetSortMode(), scorer->GetNbFold(),
                 nb_sel, ln_transf, standardize, out_path);

    featuretab_t count_tab;
    TabHeader count_tab_header(smp_info_path);
    seqVect_t feature_vect;

    EvalScore(feature_vect, count_tab_header, count_tab_path, scorer, ln_transf, standardize, idx_path);
    SortScore(feature_vect, scorer->GetSortMode());
    if (scorer->GetScoreMethod() == "ttest")
    {
        PValueAdjustmentBH(feature_vect);
    }
    ModelPrint(feature_vect, idx_path, nb_sel, count_tab_header, out_path);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
