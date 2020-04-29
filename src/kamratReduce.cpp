#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <ctime>

#include "utils/utils.hpp"
#include "utils/sample_info.hpp"
#include "reduce/utils.hpp"
#include "reduce/parse_opt_print_info.hpp"
#include "reduce/model_info.hpp"

void LoadModelInfo(std::vector<ModelInfo> &model_info_vect,
                   const std::string &kmer_count_path,
                   const SampleInfo &sample_info,
                   const size_t nb_fold,
                   evalfuncptr_t eval_func,
                   const std::string &user_method_name,
                   const std::string &transf_mode,
                   const double count_offset,
                   const size_t min_rec,
                   const size_t min_rec_abd)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        std::cerr << "ERROR: k-mer count file " << kmer_count_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;

    // Loading sample condition information //
    std::getline(kmer_count_file, line); // header line indicating sample names
    std::vector<int> header_term_nat;
    size_t sample_num = ParseHeader(header_term_nat, line, sample_info, user_method_name);
    arma::Row<size_t> sample_labels;
    ParseSampleLabel(sample_labels, header_term_nat);

    // Evaluating k-mers line by line //
    for (size_t idx_pos(kmer_count_file.tellg()); std::getline(kmer_count_file, line); idx_pos = kmer_count_file.tellg())
    {
        std::istringstream conv(line);
        double score_value;
        bool has_score_value(false);
        std::string seq, term;
        arma::mat sample_counts(1, sample_num, arma::fill::zeros);
        size_t rec(0);
        conv >> seq;
        for (size_t i_col(0), i_sample(0); conv >> term; ++i_col)
        {
            if (header_term_nat.at(i_col) >= 0)
            {
                sample_counts(0, i_sample) = static_cast<double>(std::stod(term) + count_offset);
                if (sample_counts(0, i_sample) >= min_rec_abd)
                {
                    ++rec;
                }
                if (transf_mode == "log")
                {
                    sample_counts(0, i_sample) = log(sample_counts(0, i_sample));
                }
                ++i_sample;
            }
            else if (header_term_nat.at(i_col) == -1)
            {
                score_value = std::stod(term);
                has_score_value = true;
                break;
            }
        }
        if (rec < min_rec)
        {
            continue;
        }
        if (!has_score_value)
        {
            score_value = eval_func(sample_counts, sample_labels, nb_fold, sample_info.GetNbCondition());
        }
        model_info_vect.emplace_back(seq, idx_pos, score_value);
    }

    kmer_count_file.close();
}

void ModelSort(std::vector<ModelInfo> &model_info_vect, const std::string &sort_mode)
{
    if (sort_mode == SORT_DEC_ABS)
    {
        auto comp = [](ModelInfo model1, ModelInfo model2) -> bool { return fabs(model1.GetModelScore()) > fabs(model2.GetModelScore()); };
        std::sort(model_info_vect.begin(), model_info_vect.end(), comp);
    }
    else if (sort_mode == SORT_DEC)
    {
        auto comp = [](ModelInfo model1, ModelInfo model2) -> bool { return model1.GetModelScore() > model2.GetModelScore(); };
        std::sort(model_info_vect.begin(), model_info_vect.end(), comp);
    }
    else if (sort_mode == SORT_INC_ABS)
    {
        auto comp = [](ModelInfo model1, ModelInfo model2) -> bool { return fabs(model1.GetModelScore()) < fabs(model2.GetModelScore()); };
        std::sort(model_info_vect.begin(), model_info_vect.end(), comp);
    }
    else if (sort_mode == SORT_INC)
    {
        auto comp = [](ModelInfo model1, ModelInfo model2) -> bool { return model1.GetModelScore() < model2.GetModelScore(); };
        std::sort(model_info_vect.begin(), model_info_vect.end(), comp);
    }
    else
    {
        std::cerr << "ERROR: unknown sort mode " << sort_mode << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void PValueAdjustmentBH(std::vector<ModelInfo> &model_info_vect)
{
    double tot_num = model_info_vect.size();
    for (size_t i(0); i < tot_num; ++i)
    {
        const double adjust_factor = tot_num / (i + 1);
        model_info_vect[i].ScaleScore(adjust_factor, 0, 1);
    }
    for (auto iter = model_info_vect.rbegin() + 1; iter < model_info_vect.rend(); ++iter)
    {
        double last_score = (iter - 1)->GetModelScore();
        if (iter->GetModelScore() > last_score)
        {
            iter->ScaleScore(1, last_score, last_score);
        }
    }
}

void ModelPrint(const std::vector<ModelInfo> &model_info_vect,
                const size_t top_num,
                const std::string &out_header,
                const std::string &kmer_count_path,
                const SampleInfo &sample_info)
{
    size_t parsed_top_num = (top_num == 0) ? model_info_vect.size() : top_num;
    std::ifstream kmer_count_file(kmer_count_path);
    ExitIf(!kmer_count_file.is_open(), "ERROR: cannot open count index file: " + kmer_count_path);

    std::string line;

    std::getline(kmer_count_file, line); // header line indicating sample names
    std::vector<int> header_term_nat;
    size_t sample_num = ParseHeader(header_term_nat, line, sample_info, "whatever");

    std::cout << out_header << "\t" << line.substr(line.find_first_of(" \t") + 1) << std::endl;
    for (size_t i(0); i < parsed_top_num; ++i)
    {
        std::cout << model_info_vect.at(i).GetSeq() << "\t" << model_info_vect.at(i).GetModelScore();
        GetStringLineFromDisk(line, kmer_count_file, model_info_vect.at(i).GetCountDiskPos());
        std::cout << "\t" << line.substr(line.find_first_of(" \t") + 1) << std::endl;
    }
    kmer_count_file.close();
}

int main(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string kmer_count_path, sample_info_path, eval_method(METHOD_SD), sort_mode, transf_mode;
    size_t top_num(0), min_rec(0), min_rec_abd(0);
    double count_offset(0);
    bool is_offset_set(false);
    ParseOptions(argc, argv, sample_info_path, eval_method, sort_mode, top_num, kmer_count_path, transf_mode, count_offset, min_rec, min_rec_abd, is_offset_set);
    SampleInfo sample_info;
    if (!sample_info_path.empty())
    {
        LoadSampleInfo(sample_info, sample_info_path);
    }
    if (sample_info.GetNbCondition() > 2)
    {
        std::cerr << "ERROR: sorry, the tool currently only accepts binary classification: nb_condition=" << sample_info.GetNbCondition() << std::endl;
        exit(EXIT_FAILURE);
    }

    evalfuncptr_t eval_func;
    std::vector<ModelInfo> model_info_vect;
    std::string user_method_name, out_header;
    size_t nb_fold(2);
    ParseMethod(user_method_name, eval_func, nb_fold, out_header, eval_method, sort_mode, count_offset, sample_info.GetNbCondition(), is_offset_set);
    PrintRunInfo(kmer_count_path, sample_info_path, eval_method, nb_fold, sort_mode, top_num, transf_mode, count_offset, min_rec, min_rec_abd);

    LoadModelInfo(model_info_vect, kmer_count_path, sample_info, nb_fold, eval_func, user_method_name, transf_mode, count_offset, min_rec, min_rec_abd);
    ModelSort(model_info_vect, sort_mode);
    if (eval_method == METHOD_TTEST)
    {
        PValueAdjustmentBH(model_info_vect);
    }
    ModelPrint(model_info_vect, top_num, out_header, kmer_count_path, sample_info);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
