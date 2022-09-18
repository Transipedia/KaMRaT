#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <ctime>

#include "rank_runinfo.hpp"
#include "index_loading.hpp"
#include "feature_elem.hpp"
#include "FeatureStreamer.hpp"
#include "IndexRandomAccess.hpp"
#include "scorer.hpp"


using namespace std;

using featureVect_t = std::vector<std::unique_ptr<FeatureElem>>;


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


/** Sort features idx regarding their scores.
 * @param scores A vector containing all scores. Score at position x corresponds to the xth feature.
 * @param features The features to sort.
 * @param scorer_code The scoring method to use.
 **/
void SortFeatures(const std::vector<double> & scores, std::vector<uint64_t> & features, const ScorerCode scorer_code)
{
    if (scorer_code == ScorerCode::kSNR || scorer_code == ScorerCode::kPearson || scorer_code == ScorerCode::kSpearman) // decabs
    {
        auto comp = [&scores](const uint64_t pos1, const uint64_t pos2)
            -> bool { return fabs(scores[pos1]) > fabs(scores[pos2]); };
        std::sort(features.begin(), features.end(), comp);
    }
    else if (scorer_code == ScorerCode::kTtestPi || scorer_code == ScorerCode::kDIDS || scorer_code == ScorerCode::kLR ||
             scorer_code == ScorerCode::kBayes || scorer_code == ScorerCode::kSVM ||
             scorer_code == ScorerCode::kSD || 
             scorer_code == ScorerCode::kRSD1 || scorer_code == ScorerCode::kRSD2 || scorer_code == ScorerCode::kRSD3) // dec
    {
        auto comp = [&scores](const uint64_t pos1, const uint64_t pos2)
            -> bool { return scores[pos1] > scores[pos2]; };
        std::sort(features.begin(), features.end(), comp);
    }
    else if (scorer_code == ScorerCode::kTtestPadj || scorer_code == ScorerCode::kEntropy) // inc
    {
        auto comp = [&scores](const uint64_t pos1, const uint64_t pos2)
            -> bool { return scores[pos1] < scores[pos2]; };
        std::sort(features.begin(), features.end(), comp);
    }
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


/** @param scores Sorted scores
 **/
void PrintAsIntermediate(std::vector<double> &scores, std::vector<uint64_t> &features, const size_t max_to_sel, IndexRandomAccess & ira)
{
    // --- Preprocess features/scores to free some memory ---
    uint64_t saved_space = (features.size() - max_to_sel) * sizeof(uint64_t);
    // Slice vector to remove low score features
    features.resize(max_to_sel);
    // Sort features regarding their file position
    std::sort(features.begin(), features.end());

    // --- Get matrix pointers (in place) ---
    for (uint64_t idx=0 ; idx<features.size() ; idx++) {
        features[idx] = ira.feature_to_position(features[idx]);
    }

    // --- Get file ordered features (to read the file from the beginning to the end) ---
    vector<uint64_t> feature_indexes(features.size());
    iota(std::begin(feature_indexes), std::end(feature_indexes), 0);
    auto comp = [&features](const uint64_t idx1, const uint64_t idx2)
            -> bool { return features[idx1] < features[idx2]; };
    sort(feature_indexes.begin(), feature_indexes.end(), comp);

    // --- Read from matrix file and write to stdout ---
    float * counts = new float[ira.nb_smp];
    // WARNING: Only works for constant size feature
    char * feature = new char[ira.k + 1];
    feature[ira.k] = '\0';
    for (uint64_t idx=0 ; idx<features.size() ; idx++) {
        size_t mat_idx = features[feature_indexes[idx]];
        // cout << feature_indexes[idx] << endl;
        ira.load_counts_by_file_position(mat_idx, counts, feature);
        cout << feature << "\t" << scores[idx] << "\t" << 1 << "\t";
        std::cout.write(reinterpret_cast<char *>(&mat_idx), sizeof(size_t));
        std::cout << std::endl;
    }
}


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

    IndexRandomAccess ira(idx_dir + "/idx-pos.bin", idx_dir + "/idx-mat.bin", idx_dir + "/idx-meta.bin");

    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    
    featureVect_t ft_vect;
    // after_merge = (with_path.empty() ? MakeFeatureVectFromIndex(ft_vect, idx_dir + "/idx-pos.bin", idx_mat, nb_smp, k_len)
    //                                  : MakeFeatureVectFromFile(ft_vect, with_path));
    FeatureStreamer stream = with_path.empty() ? FeatureStreamer(idx_dir + "/idx-pos.bin", idx_dir + "/idx-mat.bin", k_len, nb_smp)
                                               : FeatureStreamer(with_path);

    std::cerr << "Option parsing and metadata loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::vector<std::string> col_target_vect;
    if (rk_mthd != "sd" && rk_mthd != "rsd1" && rk_mthd != "rsd2" && rk_mthd != "rsd3" && rk_mthd != "entropy")
    {
        ParseDesign(col_target_vect, dsgn_path, colname_vect);
    }
    Scorer scorer(rk_mthd, nfold, col_target_vect);

    // Load and score all the usefull features
    size_t nb_features = 0;
    std::vector<float> count_vect;
    vector<double> scores;
    while (stream.hasNext()) {
        feature_t feature = stream.next();
        feature->EstimateCountVect(count_vect, idx_mat, nb_smp, count_mode);
        double score = scorer.EstimateScore(count_vect);
        scores.push_back(score);

        nb_features += 1;
    }

    // Postprocess variables
    after_merge = stream.merged_features;
    if (sel_top <= 0) max_to_sel = scores.size();
    else if (sel_top < 0.999999) // for avoiding when sel_top == 0.999999999999
    { max_to_sel = static_cast<size_t>(scores.size() * sel_top + 0.5); }
    else if (sel_top <= scores.size())
    { max_to_sel = static_cast<size_t>(sel_top + 0.00005); }// for avoiding when sel_top == 0.999999999999
    else
    {
        throw std::invalid_argument("number of top feature selection exceeds total feature number: " +
                                    std::to_string(static_cast<size_t>(sel_top + 0.00005)) + ">" + std::to_string(scores.size()));
    }

    // Fill a vector that will be sorted acording the scores
    std::vector<uint64_t> features(nb_features) ;
    std::iota (std::begin(features), std::end(features), 0); // Fill with 0, 1, ..., 99...
    
    // Rank the features
    SortFeatures(scores, features, scorer.GetScorerCode());

    for (int i=0 ; i<10 ; i++)
        cout << features[i] << " " << scores[features[i]] << endl;

    std::cerr << "Score evalution finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    double tot = static_cast<double>(features.size());
    if (scorer.GetScorerCode() == ScorerCode::kTtestPadj) // BH procedure
    {
        std::cerr << "\tadjusting p-values using BH procedure..." << std::endl
                  << std::endl;
        for (size_t i(features.size() - 1); i > 0; --i)
        {
            uint64_t i_score = features[i-1];
            uint64_t i1_score = features[i];

            scores[i_score] = FeatureElem::AdjustScore(
                scores[i_score],
                tot / (i + 1),
                0, scores[i1_score]
            );

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
        // TODO
        PrintWithCounts(after_merge, ft_vect, idx_mat, count_mode, nb_smp, max_to_sel);
    }
    else
    {
        // TODO
        PrintAsIntermediate(scores, features, max_to_sel, ira);
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
