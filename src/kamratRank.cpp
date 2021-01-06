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
        if (tab_header.IsColCount(i))
        {
            sum_out << tab_header.GetColNameAt(i) << "\t" << nf_vect[j++] << std::endl;
        }
    }
    sum_out.close();
}

void EvalScore(featureVect_t &feature_vect,
               std::ifstream &idx_file, const std::vector<double> nf_vect,
               std::unique_ptr<Scorer> &scorer, const TabHeader tab_header,
               const bool to_ln, const bool to_standardize, const bool no_norm)
{
    size_t nb_count = nf_vect.size();
    for (size_t i_feature(0); i_feature < feature_vect.size(); ++i_feature)
    {
        scorer->PrepareCountVect(feature_vect[i_feature], nf_vect, to_ln, to_standardize, no_norm);
        scorer->CalcFeatureStats(feature_vect[i_feature]);
        if (scorer->GetScoreMethodCode() != ScoreMethodCode::kUser) // user scoring mode but not find a score column
        {
            scorer->EvaluateScore(feature_vect[i_feature]);
        }
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
    size_t tot_num = feature_vect.size();
    for (size_t i_feature(tot_num - 2); i_feature >= 0 && i_feature < tot_num; --i_feature)
    {
        feature_vect[i_feature].AdjustScore(static_cast<double>(tot_num) / (i_feature + 1), 0, feature_vect[i_feature + 1].GetScore());
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

int main(int argc, char *argv[])
{
    std::clock_t begin_time = clock(), inter_time;
    std::string count_tab_path, smp_info_path, score_method, score_cmd, sort_mode, idx_path, nf_path, out_path;
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
    ScanCountComputeNF(feature_vect, smp_nf_vect, count_tab_header, count_tab_path, idx_path);
    PrintNF(nf_path, smp_nf_vect, count_tab_header);
    std::cerr << "Count table scanning finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("error open sample count index file: " + idx_path);
    }
    EvalScore(feature_vect, idx_file, smp_nf_vect, scorer, count_tab_header, to_ln, to_standardize);
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
