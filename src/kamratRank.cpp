#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <ctime>

#include "runinfo_files/rank_runinfo.hpp"
#include "data_struct/feature_elem.hpp"
#include "data_struct/scorer.hpp"

using featureVect_t = std::vector<std::unique_ptr<FeatureElem>>;

void LoadIndexMeta(size_t &nb_smp, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, const std::string &idx_meta_path); // in utils/index_loading.cpp
void LoadFeaturePosMap(std::unordered_map<std::string, size_t> &ft_pos_map, std::ifstream &idx_mat, const std::string &idx_pos_path,
                       const bool need_skip_code, const size_t nb_smp);                                            // in utils/index_loading.cpp
const std::string &GetTagSeq(std::string &tag_str, std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp

const bool MakeFeatureVectFromIndex(featureVect_t &ft_vect, const std::string &idx_pos_path,
                                    std::ifstream &idx_mat, const size_t nb_smp, const size_t k_len)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    size_t pos;
    std::string feature;
    while (idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t)))
    {
        if (k_len > 0)
        {
            idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t));
        }
        idx_mat.seekg(pos + nb_smp * sizeof(float)); // skip the indexed count vector
        idx_mat >> feature;
        ft_vect.emplace_back(std::make_unique<FeatureElem>(feature, pos));
    }
    ft_vect.shrink_to_fit();
    idx_pos.close();
    return false;
}

const bool MakeFeatureVectFromFile(featureVect_t &ft_vect, const std::string &with_path)
{
    bool after_merge(false), read_succ(false);
    std::ifstream with_file(with_path);
    if (!with_file.is_open())
    {
        throw std::domain_error("error open feature list file: " + with_path);
    }
    std::string feature;
    float _; // for "place-holder" column of feature scores
    size_t nb_mem_pos(1);
    std::vector<size_t> mem_pos_vect;
    std::unordered_map<std::string, size_t>::const_iterator it;
    while (with_file >> feature >> _ >> nb_mem_pos)
    {
        read_succ = true;
        mem_pos_vect.resize(nb_mem_pos);
        with_file.ignore(1);
        with_file.read(reinterpret_cast<char *>(&mem_pos_vect[0]), nb_mem_pos * sizeof(size_t));
        if (nb_mem_pos > 0)
        {
            after_merge = true;
        }
        ft_vect.emplace_back(std::make_unique<FeatureElem>(feature, mem_pos_vect));
    }
    ft_vect.shrink_to_fit();
    with_file.close();
    if (!read_succ)
    {
        throw std::invalid_argument("not valid file for input sequences: " + with_path);
    }
    return after_merge;
}

void ParseDesign(std::vector<std::string> &col_target_vect, const std::string &dsgn_path, const std::vector<std::string> &colname_vect)
{
    const size_t nb_smp = colname_vect.size() - 1; // set aside the first column representing features
    
    std::ifstream dsgn_file(dsgn_path);
    if (!dsgn_file.is_open())
    {
        throw std::invalid_argument("error open design file: " + dsgn_path);
    }
    std::map<std::string, size_t> col_name2num;
    for (size_t i(1); i <= nb_smp; ++i)
    {
        if (!col_name2num.insert({colname_vect[i], i - 1}).second)
        {
            throw std::invalid_argument("duplicate column name in the matrix: " + colname_vect[i]);
        }
    }

    col_target_vect.resize(nb_smp, "");
    std::string line, smp_name, condi("");
    std::istringstream line_conv;
    while (std::getline(dsgn_file, line))
    {
        line_conv.str(line);
        if (!(line_conv >> smp_name >> condi))
        {
            throw std::domain_error("failed in design file parsing: line " + line);
        }
        line_conv.clear();
        auto iter = col_name2num.find(smp_name);
        if (iter == col_name2num.cend())
        {
            std::cerr << "[info] a sample in the design file not found in the matrix header line: " + smp_name << std::endl;
            continue;
        }
        col_target_vect[iter->second] = condi;
    }
    for (size_t i(0); i < nb_smp; ++i)
    {
        if (col_target_vect[i].empty())
        {
            throw std::invalid_argument("a column in the matrix was not well annotated by the design file: " + colname_vect[i + 1]);
        }
        // std::cerr << colname_vect[i + 1] << "\t" << col_target_vect[i] << std::endl; // for check
    }
    dsgn_file.close();
}

void SortScore(featureVect_t &ft_vect, const ScorerCode scorer_code)
{
    if (scorer_code == ScorerCode::kSNR || scorer_code == ScorerCode::kPearson || scorer_code == ScorerCode::kSpearman) // decabs
    {
        auto comp = [](const std::unique_ptr<FeatureElem> &ft1, const std::unique_ptr<FeatureElem> &ft2)
            -> bool { return fabs(ft1->GetScore()) > fabs(ft2->GetScore()); };
        std::sort(ft_vect.begin(), ft_vect.end(), comp);
    }
    else if (scorer_code == ScorerCode::kTtestPi || scorer_code == ScorerCode::kDIDS || scorer_code == ScorerCode::kLR ||
             scorer_code == ScorerCode::kBayes || scorer_code == ScorerCode::kSVM ||
             scorer_code == ScorerCode::kSD || 
             scorer_code == ScorerCode::kRSD1 || scorer_code == ScorerCode::kRSD2 || scorer_code == ScorerCode::kRSD3) // dec
    {
        auto comp = [](const std::unique_ptr<FeatureElem> &ft1, const std::unique_ptr<FeatureElem> &ft2)
            -> bool { return ft1->GetScore() > ft2->GetScore(); };
        std::sort(ft_vect.begin(), ft_vect.end(), comp);
    }
    else if (scorer_code == ScorerCode::kTtestPadj || scorer_code == ScorerCode::kEntropy) // inc
    {
        auto comp = [](const std::unique_ptr<FeatureElem> &ft1, const std::unique_ptr<FeatureElem> &ft2)
            -> bool { return ft1->GetScore() < ft2->GetScore(); };
        std::sort(ft_vect.begin(), ft_vect.end(), comp);
    }
    // else if (...) // incabs, useless
    // {
    //     auto comp = [](const std::unique_ptr<FeatureElem> &ft1, const std::unique_ptr<FeatureElem> &ft2)
    //         -> bool { return fabs(ft1->GetScore()) < fabs(ft2->GetScore()); };
    //     std::sort(ft_vect.begin(), ft_vect.end(), comp);
    // }
}

void PrintHeader(const bool after_merge, const std::vector<std::string> &colname_vect, const std::string &scorer_name)
{
    if (after_merge)
    {
        std::cout << "contig\tnb-merged-kmer"
                  << "\t";
    }
    std::cout << colname_vect[0] << "\t" << scorer_name;
    for (size_t i_col(1); i_col < colname_vect.size(); ++i_col)
    {
        std::cout << "\t" << colname_vect[i_col];
    }
    std::cout << std::endl;
}

void PrintWithCounts(const bool after_merge, const featureVect_t &ft_vect, std::ifstream &idx_mat,
                     const std::string &count_mode, const size_t nb_smp, const size_t max_to_sel)
{
    std::vector<float> count_vect;
    std::string rep_seq;
    for (size_t i_ft(0); i_ft < max_to_sel; ++i_ft)
    {
        std::cout << ft_vect[i_ft]->GetFeature();
        if (after_merge)
        {
            std::cout << "\t" << ft_vect[i_ft]->GetNbMemPos();
            std::cout << "\t" << GetTagSeq(rep_seq, idx_mat, ft_vect[i_ft]->GetRepPos(), nb_smp);
        }
        std::cout << "\t" << ft_vect[i_ft]->GetScore();
        ft_vect[i_ft]->EstimateCountVect(count_vect, idx_mat, nb_smp, count_mode);
        for (float x : count_vect)
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
    }
}

void PrintAsIntermediate(const featureVect_t &ft_vect, const size_t max_to_sel)
{
    size_t p;
    for (size_t i_ft(0); i_ft < max_to_sel; ++i_ft)
    {
        std::cout << ft_vect[i_ft]->GetFeature() << "\t" << ft_vect[i_ft]->GetScore() << "\t"
                  << ft_vect[i_ft]->GetNbMemPos() << "\t";
        p = ft_vect[i_ft]->GetRepPos();
        std::cout.write(reinterpret_cast<char *>(&p), sizeof(size_t));
        std::cout << std::endl;
    }
}

using namespace std;

int RankMain(int argc, char *argv[])
{
    RankWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, rk_mthd, with_path, count_mode("rep"), dsgn_path, out_path;
    float sel_top(-1); // negative value means without selection, print all features
    size_t nfold, nb_smp, k_len, max_to_sel;
    bool with_counts(false), after_merge(false), _stranded; // _stranded not needed in KaMRaT-rank
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, rk_mthd, nfold, with_path, count_mode, dsgn_path, sel_top, out_path, with_counts);
    PrintRunInfo(idx_dir, rk_mthd, nfold, with_path, count_mode, dsgn_path, sel_top, out_path, with_counts);
    LoadIndexMeta(nb_smp, k_len, _stranded, colname_vect, idx_dir + "/idx-meta.bin");

    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    featureVect_t ft_vect;
    after_merge = (with_path.empty() ? MakeFeatureVectFromIndex(ft_vect, idx_dir + "/idx-pos.bin", idx_mat, nb_smp, k_len)
                                     : MakeFeatureVectFromFile(ft_vect, with_path));
    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    if (sel_top <= 0)
    {
        max_to_sel = ft_vect.size();
    }
    else if (sel_top < 0.999999) // for avoiding when sel_top == 0.999999999999
    {
        max_to_sel = static_cast<size_t>(ft_vect.size() * sel_top + 0.5);
    }
    else if (sel_top <= ft_vect.size())
    {
        max_to_sel = static_cast<size_t>(sel_top + 0.00005); // for avoiding when sel_top == 0.999999999999
    }
    else
    {
        throw std::invalid_argument("number of top feature selection exceeds total feature number: " +
                                    std::to_string(static_cast<size_t>(sel_top + 0.00005)) + ">" + std::to_string(ft_vect.size()));
    }
    std::vector<std::string> col_target_vect;
    if (rk_mthd != "sd" && rk_mthd != "rsd1" && rk_mthd != "rsd2" && rk_mthd != "rsd3" && rk_mthd != "entropy")
    {
        ParseDesign(col_target_vect, dsgn_path, colname_vect);
    }

    // --- Score the vectors ---
    cout << "--------------------------------------------------------------" << endl;
    Scorer scorer(rk_mthd, nfold, col_target_vect);
    cout << "ft vect size " << ft_vect.size() << endl;
    
    std::vector<std::vector<float> > count_vects; // Used to load multiple features
    uint ftx_step = 100; // Number of feature to read at onece
    for (uint ftx=0 ; ftx<ft_vect.size() ; ftx += ftx_step) {
        // 1 - Loading data (sequencial)
        for (uint offset=0 ; offset<ftx_step ; offset++) {
            ft_vect[ftx + offset]->EstimateCountVect(count_vects[offset], idx_mat, nb_smp, count_mode);    
        }
        // 2 - Estimate score (parallel)
        std::vector<double> scores = scorer.EstimateScores(count_vects);
        // 3 - Set the score (sequencial)
        for (uint offset=0 ; offset<ftx_step ; offset++) {
            ft_vect[ftx + offset]->SetScore(scores[offset]);
        }
    }
    cout << "//////////////////////////////////////////////////////////////" << endl;

    // std::vector<float> count_vect;
    // for (const auto &ft : ft_vect)
    // {
    //     ft->EstimateCountVect(count_vect, idx_mat, nb_smp, count_mode);
    //     ft->SetScore(scorer.EstimateScore(count_vect));
    // }
    SortScore(ft_vect, scorer.GetScorerCode());
    std::cerr << "Score evalution finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();
    if (scorer.GetScorerCode() == ScorerCode::kTtestPadj) // BH procedure
    {
        std::cerr << "\tadjusting p-values using BH procedure..." << std::endl
                  << std::endl;
        for (size_t tot_num(ft_vect.size()), i_ft(tot_num - 2); i_ft >= 0 && i_ft < tot_num; --i_ft)
        {
            ft_vect[i_ft]->AdjustScore(static_cast<double>(tot_num) / (i_ft + 1), 0, ft_vect[i_ft + 1]->GetScore());
        }
        std::cerr << "P-value adjusting finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
        inter_time = clock();
    }

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
    if (with_counts)
    {
        PrintHeader(after_merge, colname_vect, scorer.GetScorerName());
        PrintWithCounts(after_merge, ft_vect, idx_mat, count_mode, nb_smp, max_to_sel);
    }
    else
    {
        PrintAsIntermediate(ft_vect, max_to_sel);
    }
    idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << max_to_sel << " features have been written." << std::endl;
    std::cerr << "Output finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
