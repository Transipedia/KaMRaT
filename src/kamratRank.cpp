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

#include "common/tab_header.hpp"
#include "rank/rank_runinfo.hpp" // including scorer, scorer includes feature_elem

const void ScanCountComputeNF(featureVect_t &feature_vect, std::vector<double> &nf_vect, TabHeader &tab_header,
                              const std::string &raw_counts_path, const std::string &idx_path, const std::string &rep_column)
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

    std::string line;
    std::getline(kmer_count_instream, line);
    std::istringstream conv(line);
    tab_header.MakeColumns(conv, rep_column);
    conv.clear();

    nf_vect.resize(tab_header.GetNbCount(), 0);
    std::cerr << "\t => Number of sample parsed: " << nf_vect.size() << std::endl;

    std::ofstream idx_file(idx_path);
    if (!idx_file.is_open()) // to ensure the file is opened
    {
        throw std::domain_error("error open file: " + idx_path);
    }

    std::vector<float> count_vect;
    std::string value_str;
    double mean_sample_sum(0), rep_value;
    while (std::getline(kmer_count_instream, line))
    {
        conv.str(line);
        rep_value = tab_header.ParseRowStr(count_vect, value_str, conv);
        feature_vect.emplace_back(rep_value, count_vect, value_str, idx_file);
        for (size_t i(0); i < count_vect.size(); ++i)
        {
            nf_vect[i] += count_vect[i];
            mean_sample_sum += count_vect[i];
        }
        count_vect.clear();
        conv.clear();
    }
    mean_sample_sum /= nf_vect.size(); // mean of sample-sum vector
    for (size_t i(0); i < nf_vect.size(); ++i)
    {
        if (nf_vect[i] == 0)
        {
            nf_vect[i] = 0;
        }
        else
        {
            nf_vect[i] = mean_sample_sum / nf_vect[i];
        }
    }
    raw_counts_file.close();
    idx_file.close();
}

void PrintNF(const std::string &smp_sum_outpath, const std::vector<double> &nf_vect, const TabHeader &tab_header)
{
    std::ofstream sum_out(smp_sum_outpath);
    for (size_t i(1), j(0); i < tab_header.GetNbCol(); ++i)
    {
        if (tab_header.IsColCount(i))
        {
            sum_out << tab_header.GetColNameAt(i) << "\t" << nf_vect[j++] << std::endl;
        }
    }
    sum_out.close();
}

void EvalScore(featureVect_t &feature_vect,
               std::ifstream &idx_file, const std::vector<double> &nf_vect,
               std::unique_ptr<Scorer> &scorer, const TabHeader &tab_header,
               const bool to_ln, const bool to_standardize, const bool no_norm)
{
    std::vector<size_t> label_vect;
    scorer->LoadSampleLabels(tab_header.GetSampleLabelVect(label_vect));
    for (size_t i_feature(0); i_feature < feature_vect.size(); ++i_feature)
    {
        scorer->PrepareCountVect(feature_vect[i_feature], nf_vect, idx_file, to_ln, to_standardize, no_norm);
        scorer->CalcFeatureStats(feature_vect[i_feature]);
        feature_vect[i_feature].SetScore(scorer->EvaluateScore(feature_vect[i_feature]));
    }
}

void SortScore(featureVect_t &feature_vect, const SortModeCode sort_mode_code)
{
    if (sort_mode_code == SortModeCode::kDecAbs)
    {
        auto comp = [](const FeatureElem &feature1, const FeatureElem &feature2) -> bool { return fabs(feature1.GetScore()) > fabs(feature2.GetScore()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode_code == SortModeCode::kDec)
    {
        auto comp = [](const FeatureElem &feature1, const FeatureElem &feature2) -> bool { return feature1.GetScore() > feature2.GetScore(); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode_code == SortModeCode::kIncAbs)
    {
        auto comp = [](const FeatureElem &feature1, const FeatureElem &feature2) -> bool { return fabs(feature1.GetScore()) < fabs(feature2.GetScore()); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
    else if (sort_mode_code == SortModeCode::kInc)
    {
        auto comp = [](const FeatureElem &feature1, const FeatureElem &feature2) -> bool { return feature1.GetScore() < feature2.GetScore(); };
        std::sort(feature_vect.begin(), feature_vect.end(), comp);
    }
}

void PValueAdjustmentBH(featureVect_t &feature_vect)
{
    size_t tot_num = feature_vect.size();
    for (size_t i_feature(tot_num - 2); i_feature >= 0 && i_feature < tot_num; --i_feature)
    {
        feature_vect[i_feature].AdjustScore(static_cast<double>(tot_num) / (i_feature + 1), 0, feature_vect[i_feature + 1].GetScore());
    }
}

void PrintHeader(std::ostream &out_s, const TabHeader &tab_header, const ScoreMethodCode score_method_code)
{
    std::cout << tab_header.GetColNameAt(0);
    if (score_method_code == ScoreMethodCode::kRelatSD) // relative sd ranking output stats for all samples
    {
        std::cout << "\t" << kScoreMethodName[score_method_code] << "\tmean\tsd";
    }
    else if (score_method_code != ScoreMethodCode::kUser) // user ranking do not output score or condition stats
    {
        std::cout << "\t" << kScoreMethodName[score_method_code];
        for (size_t i(0); i < tab_header.GetNbCondition(); ++i)
        {
            std::cout << "\tmean" << static_cast<char>('A' + i);
        }
        for (size_t i(0); i < tab_header.GetNbCondition(); ++i)
        {
            std::cout << "\tsd" << static_cast<char>('A' + i);
        }
    }
    std::string value_str;
    for (size_t i(1); i < tab_header.GetNbCol(); ++i)
    {
        if (!tab_header.IsColCount(i)) // first output non-sample values
        {
            out_s << "\t" << tab_header.GetColNameAt(i);
        }
        else
        {
            value_str += ("\t" + tab_header.GetColNameAt(i));
        }
    }
    out_s << value_str << std::endl; // then output sample counts
}

void PrintFeature(std::ostream &out_s, const FeatureElem &feature_elem, std::ifstream &idx_file,
                  const size_t nb_count, const ScoreMethodCode score_method_code)
{
    static std::vector<float> count_vect;
    static std::string value_str;
    feature_elem.RetrieveCountVect(count_vect, idx_file, nb_count);
    feature_elem.RetrieveValueStr(value_str, idx_file, nb_count);
    size_t split_pos = value_str.find_first_of(" \t");
    out_s << value_str.substr(0, split_pos);
    if (score_method_code != ScoreMethodCode::kUser)
    {
        std::cout << "\t" << feature_elem.GetScore();
        for (const double m : feature_elem.GetCondiMeanVect())
        {
            std::cout << "\t" << m;
        }
        for (const double s : feature_elem.GetCondiStddevVect())
        {
            std::cout << "\t" << s;
        }
    }
    if (split_pos != std::string::npos) // if some other value remains in value string
    {
        out_s << value_str.substr(split_pos);
    }
    for (float c : count_vect)
    {
        out_s << "\t" << c;
    }
    out_s << std::endl;
}

void ModelPrint(featureVect_t &feature_vect, std::ifstream &idx_file, const size_t nb_sel,
                const TabHeader &tab_header, const std::string &out_path, const ScoreMethodCode score_method_code)
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
    PrintHeader(std::cout, tab_header, score_method_code);
    size_t parsed_nb_sel = (nb_sel == 0 ? feature_vect.size() : nb_sel);
    for (size_t i(0); i < parsed_nb_sel; ++i)
    {
        PrintFeature(std::cout, feature_vect[i], idx_file, tab_header.GetNbCount(), score_method_code);
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
    std::string count_tab_path, smp_info_path, idx_path, nf_path, out_path;
    bool to_ln(false), to_standardize(false), no_norm(false);
    size_t nb_sel(0);
    std::unique_ptr<Scorer> scorer;

    ParseOptions(argc, argv, idx_path, nf_path, smp_info_path, scorer, nb_sel, to_ln, to_standardize, no_norm, out_path, count_tab_path);
    PrintRunInfo(count_tab_path, idx_path, nf_path, smp_info_path, scorer, nb_sel, to_ln, to_standardize, no_norm, out_path);
    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    TabHeader count_tab_header(smp_info_path);
    featureVect_t feature_vect;
    std::vector<double> smp_nf_vect;
    if (scorer->GetScoreMethodCode() != ScoreMethodCode::kUser)
    {
        ScanCountComputeNF(feature_vect, smp_nf_vect, count_tab_header, count_tab_path, idx_path, "");
    }
    else
    {
        ScanCountComputeNF(feature_vect, smp_nf_vect, count_tab_header, count_tab_path, idx_path, scorer->GetRepColname());
    }
    std::cerr << "Count table scanning finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();
    if (!no_norm)
    {
        PrintNF(nf_path, smp_nf_vect, count_tab_header);
    }

    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("error open sample count index file: " + idx_path);
    }
    if (scorer->GetScoreMethodCode() != ScoreMethodCode::kUser)
    {
        EvalScore(feature_vect, idx_file, smp_nf_vect, scorer, count_tab_header, to_ln, to_standardize, no_norm);
    }
    std::cerr << "Score evalution finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();
    SortScore(feature_vect, scorer->GetSortModeCode());
    std::cerr << "Score sorting finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();
    if (scorer->GetScoreMethodCode() == ScoreMethodCode::kTtest)
    {
        PValueAdjustmentBH(feature_vect);
        std::cerr << "p-value adjusted by BH method..." << std::endl
                  << std::endl;
        std::cerr << "P-value adjusting finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
        inter_time = clock();
    }
    ModelPrint(feature_vect, idx_file, nb_sel, count_tab_header, out_path, scorer->GetScoreMethodCode());
    idx_file.close();

    std::cerr << "Output finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;

    return EXIT_SUCCESS;
}
