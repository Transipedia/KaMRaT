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


/**
  * @param scores The score of each feature according to the input order
  * @param features The index of the features sorted with the best score first
  * @param max_to_sel Maximum number of element to keep in the output
  * @param stream Object that will allow to enumerate features onr by one in the input order
  * @param idx_mat matrix indexes file
  * @param nb_smp number of columns in the matrix
  * @param count_mode Counting mode
 **/
void PrintWithCounts_features(const std::vector<double> & scores, std::vector<uint64_t> & features, size_t max_to_sel, FeatureStreamer & stream, ifstream & idx_mat, size_t nb_smp, std::string count_mode)
{
    // for (double val : scores)
    //     std::cout << val << " ";
    // std::cout << std::endl;
    std::cerr << "max to sel " << max_to_sel << std::endl;
    // Resize the feature vector to only keep the max_to_sel best ones
    features.resize(max_to_sel);
    std::cerr << "features size " << features.size() << endl;
    for (uint64_t feature : features)
        std::cerr << feature << " ";
    std::cerr << std::endl;
    std::sort(features.begin(), features.end());
    std::cerr << "features size " << features.size() << endl;
    for (uint64_t feature : features)
        std::cerr << feature << " ";
    std::cerr << std::endl;

    std::vector<float> count_vect;

    uint64_t feature_idx = 0, idx = 0;
    while (stream.hasNext() and feature_idx < max_to_sel) {
        feature_t feature = stream.next();

        // If the next feature of interest is not yet reached
        if (features[feature_idx] > idx++) {
            continue;
        }
        std::cerr << "In " << features[feature_idx] << std::endl;
        feature->EstimateCountVect(count_vect, idx_mat, nb_smp, count_mode);

        // Print the current feature
        std::string rep_seq;
        std::cout << feature->GetFeature() << "\t" << feature->GetNbMemPos();
        std::cout << "\t" << GetTagSeq(rep_seq, idx_mat, feature->GetRepPos(), nb_smp);
        std::cout << "\t" << scores[features[feature_idx]];
        for (float x : count_vect)
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
    
        feature_idx += 1;
    }
}


void PrintWithCounts_kmers(const std::vector<double> &scores, const std::vector<uint64_t> &features,
                     std::ifstream &idx_mat, const std::string &count_mode, const size_t nb_smp,
                     const size_t max_to_sel, IndexRandomAccess & ira)
{
    float * counts = new float[ira.nb_smp];
    char * feature = new char[ira.k + 1];
    feature[ira.k] = '\0';

    for (size_t idx(0); idx < max_to_sel; ++idx)
    {
        // WARNING: Only works for kmer features
        size_t mat_idx = ira.feature_to_position(features[idx]);
        ira.load_counts_by_file_position(mat_idx, counts, feature);

        std::cout << feature;
        std::cout << "\t" << scores[features[idx]];
        
        for (size_t idx(0) ; idx<ira.nb_smp ; idx++)
        {
            std::cout << "\t" << counts[idx];
        }
        std::cout << std::endl;
    }

    delete[] counts;
    delete[] feature;
}


typedef struct feature_pos_s {
    uint64_t feature;
    uint64_t file_pos;
} feature_pos_t;



void PrintAsIntermediate_features(std::vector<double> & scores, std::vector<uint64_t> & features, size_t max_to_sel, FeatureStreamer & stream, ifstream & idx_mat, size_t nb_smp, std::string count_mode)
{
    features.resize(max_to_sel);
    std::sort(features.begin(), features.end());

    std::vector<float> count_vect;

    uint64_t feature_idx = 0, idx = 0;
    while (stream.hasNext()) {
        feature_t feature = stream.next();

        // If the next feature of interest is not yet reached
        if (features[feature_idx] > idx++) {
            continue;
        }
        feature->EstimateCountVect(count_vect, idx_mat, nb_smp, count_mode);

        // Print the current feature
        std::cout << feature->GetFeature() << "\t" << scores[features[feature_idx]] << "\t"
                  << feature->GetNbMemPos() << "\t";
        size_t p = feature->GetRepPos();
        std::cout.write(reinterpret_cast<char *>(&p), sizeof(size_t));
        std::cout << std::endl;

        feature_idx += 1;
    }
}

/** Print the outputs. Features need to be reloaded from the matrix. This function is highly
  * modified to read all the files in the right order (begin to end and not random access) 
  * @param scores Scores by feature. scores[i] is the score for the ith feature.
  * @param features Feature idxs sorted by score.
  * @param max_to_sel Number max of features to keep in the output.
  * @param ira Object to allow index random accesses.
 **/
void PrintAsIntermediate_kmers(std::vector<double> &scores, std::vector<uint64_t> &features, const size_t max_to_sel, IndexRandomAccess & ira)
{
    // --- Preprocess features/scores to free some memory ---
    uint64_t saved_space = (features.size() - max_to_sel) * sizeof(uint64_t);
    // Slice vector to remove low score features
    features.resize(max_to_sel);
    // Sort features regarding their file position
    std::sort(features.begin(), features.end());

    // --- Get matrix pointers ---
    std::vector<feature_pos_t> feature_positions(max_to_sel);
    for (uint64_t idx=0 ; idx<features.size() ; idx++) {
        feature_positions[idx] = {features[idx], ira.feature_to_position(features[idx])};
    }

    // --- Get file ordered features (to read the file from the beginning to the end) ---
    auto comp = [](const feature_pos_t f1, const feature_pos_t f2)
            -> bool { return f1.file_pos < f2.file_pos; };
    sort(feature_positions.begin(), feature_positions.end(), comp);

    // --- Read from matrix file and write to stdout ---
    float * counts = new float[ira.nb_smp];
    char * feature = new char[ira.k + 1];
    feature[ira.k] = '\0';
    for (uint64_t idx=0 ; idx<feature_positions.size() ; idx++) {
        size_t mat_idx = feature_positions[idx].file_pos;
        ira.load_counts_by_file_position(mat_idx, counts, feature);
        std::cout << feature << "\t" << scores[feature_positions[idx].feature] << "\t" << 1 << "\t";
        std::cout.write(reinterpret_cast<char *>(&mat_idx), sizeof(size_t));
        std::cout << std::endl;
    }
    
    delete[] counts;
    delete[] feature;
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
        if (after_merge) {
            std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
            FeatureStreamer stream(with_path);
            PrintWithCounts_features(scores, features, max_to_sel, stream, idx_mat, nb_smp, count_mode);
        } else {
            PrintWithCounts_kmers(scores, features, idx_mat, count_mode, nb_smp, max_to_sel, ira);
        }
    }
    else
    {
        if (with_path.empty())
            PrintAsIntermediate_kmers(scores, features, max_to_sel, ira);
        else {
            std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
            FeatureStreamer stream(with_path);
            PrintAsIntermediate_features(scores, features, max_to_sel, stream, idx_mat, nb_smp, count_mode);
        }
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
