#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <cassert>
#include <ctime>

#include <getopt.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <mlpack/core.hpp>

#include "utils.h"
#include "kmer_model.h"
#include "eval_methods.h"

const size_t MetaDataParser(std::unordered_map<std::string, size_t> &sample_info,
                            const std::string &sample_info_path)
{
    std::ifstream sample_condition_file(sample_info_path);
    if (!sample_condition_file.is_open())
    {
        std::cerr << "ERROR: meta data file " << sample_info_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::unordered_map<std::string, size_t> condition_dict;
    size_t class_num(0), label;
    while (std::getline(sample_condition_file, line))
    {
        std::istringstream conv(line);
        std::string sample, condition;
        conv >> sample >> condition;
        if (conv.fail()) // dealing with the case when sample_info contains only sample column
        {
            condition = "condition is abscent";
        }

        auto iter = condition_dict.find(condition);
        if (iter == condition_dict.cend())
        {
            condition_dict.insert({condition, class_num});
            label = class_num;
            class_num++;
        }
        else
        {
            label = iter->second;
        }

        if (!sample_info.insert({sample, label}).second)
        {
            std::cerr << "ERROR: meta data file has duplicated sample name or sample-condition pair: "
                      << "(" << sample << ", " << condition << ")." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    sample_condition_file.close();

    if (class_num > 2)
    {
        std::cerr << "ERROR: sorry, the tool currently only accepts binary classification..." << std::endl;
        exit(EXIT_FAILURE);
    }

    return class_num;
}

void MethodParser(std::string &user_method_name,
                  evalfuncptr_t &eval_func,
                  size_t &fold_num,
                  std::string &out_header,
                  std::string &eval_method,
                  std::string &sort_mode,
                  const size_t class_num)
{
    size_t split_pos = eval_method.find(":");
    std::string parsed_eval_method = eval_method.substr(0, split_pos);

    out_header = "feature\t";
    if (parsed_eval_method == METHOD_NB)
    {
        if (class_num < 2)
        {
            std::cerr << "ERROR: naive Bayes evaluation cannot be applied on single-class data." << std::endl;
            exit(EXIT_FAILURE);
        }
        eval_func = nb_evaluation;
        if (split_pos != std::string::npos)
        {
            fold_num = std::stoi(eval_method.substr(split_pos + 1));
        }
        out_header += "f1.naiveBayes";
        eval_method = parsed_eval_method;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (parsed_eval_method == METHOD_LR)
    {
        if (class_num < 2)
        {
            std::cerr << "ERROR: logistic regression evaluation cannot be applied on single-class data." << std::endl;
            exit(EXIT_FAILURE);
        }
        eval_func = lr_evaluation;
        if (split_pos != std::string::npos)
        {
            fold_num = std::stoi(eval_method.substr(split_pos + 1));
        }
        out_header += "f1.LogitReg";
        eval_method = parsed_eval_method;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (parsed_eval_method == METHOD_SD)
    {
        eval_func = sd_evaluation;
        out_header += "sd";
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (parsed_eval_method == METHOD_RSD)
    {
        eval_func = rsd_evaluation;
        out_header += "relat_sd";
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (parsed_eval_method == METHOD_SDC)
    {
        if (class_num < 2)
        {
            std::cerr << "ERROR: absolute contrast of standard deviation evaluation cannot be applied on single-class data." << std::endl;
            exit(EXIT_FAILURE);
        }
        eval_func = sdc_evaluation;
        out_header += "sd.abs.ctrst";
        if (sort_mode.empty())
        {
            sort_mode = std::string(SORT_DEC) + ':' + std::string(SORT_ABS);
        }
    }
    else if (parsed_eval_method == METHOD_TTEST)
    {
        if (class_num < 2)
        {
            std::cerr << "ERROR: t-test evaluation cannot be applied on single-class data." << std::endl;
            exit(EXIT_FAILURE);
        }
        eval_func = ttest_evaluation;
        out_header += "pval.t-test";
        if (sort_mode.empty())
        {
            sort_mode = SORT_INC;
        }
    }
    else if (parsed_eval_method == METHOD_ES)
    {
        eval_func = es_evaluation;
        out_header += "effect_size";
        if (sort_mode.empty())
        {
            sort_mode = std::string(SORT_DEC) + ':' + std::string(SORT_ABS);
        }
    }
    else if (parsed_eval_method == METHOD_USER)
    {
        if (split_pos == std::string::npos)
        {
            std::cerr << "ERROR: missing user-defined method name." << std::endl
                      << "       please type kamratReduce -h for help..." << std::endl;
            exit(EXIT_FAILURE);
        }
        user_method_name = eval_method.substr(split_pos + 1);
        eval_func = NULL;
        out_header += user_method_name;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else
    {
        std::cerr << "ERROR: unknown evaluation method name (" << eval_method << ")." << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

const size_t SampleLabelParser(arma::Row<size_t> &sample_labels,
                               const std::string &header_line,
                               const std::unordered_map<std::string, size_t> &sample_info)
{
    std::istringstream conv(header_line);
    std::string sample_x;
    std::vector<size_t> counts;

    conv >> sample_x;
    for (conv >> sample_x; !conv.fail(); conv >> sample_x)
    {
        auto iter = sample_info.find(sample_x); // to check if the current column is a sample column
        if (iter != sample_info.cend())
        {
            counts.push_back(iter->second);
        }
        else if (sample_info.empty())
        {
            counts.push_back(0);
        }
    }

    size_t n_sample = counts.size();
    sample_labels.set_size(n_sample);
    sample_labels = arma::conv_to<arma::Row<size_t>>::from(counts);

    return n_sample;
}

void KmerEvaluate(std::vector<KmerModel> &kmer_model_vec,
                  const std::string &kmer_count_path,
                  const std::unordered_map<std::string, size_t> &sample_info,
                  const size_t fold_num,
                  evalfuncptr_t eval_func,
                  const std::string &user_method_name)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        std::cerr << "ERROR: k-mer count file " << kmer_count_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t pos = kmer_count_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && kmer_count_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    // Loading sample condition information //
    arma::Row<size_t> sample_labels;
    std::string header_line;
    std::getline(kmer_count_instream, header_line); // header line indicating sample names
    size_t sample_num = SampleLabelParser(sample_labels, header_line, sample_info);

    // Evaluating k-mers line by line //
    std::string count_line;
    size_t n_line(0);
    while (std::getline(kmer_count_instream, count_line))
    {
        if (n_line > 0 && n_line % 10000000 == 0)
        {
            std::cerr << "\t" << n_line << " k-mers has been processed..." << std::endl;
        }
        KmerModel kmer_model(sample_num, count_line, header_line, sample_labels, sample_info, fold_num, eval_func, user_method_name);
        kmer_model_vec.push_back(kmer_model);
        n_line++;
    }

    if (kmer_count_file.is_open())
    {
        kmer_count_file.close();
    }
}

void KmerSort(std::vector<KmerModel> &kmer_model_list, const std::string &sort_mode)
{
    size_t split_pos = sort_mode.find(":");
    std::string order, sub_func;

    if (split_pos != std::string::npos)
    {
        order = sort_mode.substr(0, split_pos);
        sub_func = sort_mode.substr(split_pos + 1);
    }
    else
    {
        order = sort_mode;
    }

    if (!sub_func.empty() && sub_func != SORT_ABS)
    {
        std::cerr << "ERROR: unknown sort sub-function " << sort_mode << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (order == SORT_DEC && sub_func == SORT_ABS)
    {
        auto comp = [](KmerModel kmer1, KmerModel kmer2) -> bool { return fabs(kmer1.GetScore()) > fabs(kmer2.GetScore()); };
        std::sort(kmer_model_list.begin(), kmer_model_list.end(), comp);
    }
    else if (order == SORT_DEC)
    {
        auto comp = [](KmerModel kmer1, KmerModel kmer2) -> bool { return kmer1.GetScore() > kmer2.GetScore(); };
        std::sort(kmer_model_list.begin(), kmer_model_list.end(), comp);
    }
    else if (order == SORT_INC && sub_func == SORT_ABS)
    {
        auto comp = [](KmerModel kmer1, KmerModel kmer2) -> bool { return fabs(kmer1.GetScore()) < fabs(kmer2.GetScore()); };
        std::sort(kmer_model_list.begin(), kmer_model_list.end(), comp);
    }
    else if (order == SORT_INC)
    {
        auto comp = [](KmerModel kmer1, KmerModel kmer2) -> bool { return kmer1.GetScore() < kmer2.GetScore(); };
        std::sort(kmer_model_list.begin(), kmer_model_list.end(), comp);
    }
    else
    {
        std::cerr << "ERROR: unknown sort mode " << sort_mode << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void PValue_Adjustment_BH(std::vector<KmerModel> &kmer_model_list)
{
    double tot_num = kmer_model_list.size();

    for (size_t i = 0; i < tot_num; i++)
    {
        const double adjust_factor = tot_num / (i + 1);
        kmer_model_list[i].ScaleScore(adjust_factor, 0, 1);
    }
    for (auto iter = kmer_model_list.rbegin() + 1; iter < kmer_model_list.rend(); iter++)
    {
        double last_score = (iter - 1)->GetScore();
        if (iter->GetScore() > last_score)
        {
            iter->ScaleScore(1, last_score, last_score);
        }
    }
}

void KmerPrint(const std::vector<KmerModel> &kmer_model_list, const size_t top_num, const std::string &out_header)
{
    size_t parsed_top_num;

    if (top_num == 0)
    {
        parsed_top_num = kmer_model_list.size();
    }
    else
    {
        parsed_top_num = top_num;
    }

    std::cout << out_header << std::endl;
    for (size_t i = 0; i < parsed_top_num; i++)
    {
        kmer_model_list[i].Print();
    }
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock();
    std::string kmer_count_path, sample_info_path, eval_method(METHOD_NB), user_method_name, sort_mode, out_header;
    size_t fold_num(2), top_num(0), class_num(1);
    bool select_required(false), binary_input(false), has_condition(true);
    evalfuncptr_t eval_func;
    std::vector<KmerModel> kmer_model_vec;
    std::unordered_map<std::string, size_t> sample_info;

    int opt;
    while ((opt = getopt(argc, argv, "d:m:s:n:h")) != -1)
    {
        switch (opt)
        {
        case 'd':
            sample_info_path = optarg;
            break;
        case 'm':
            eval_method = optarg;
            break;
        case 's':
            sort_mode = optarg;
            break;
        case 'n':
            top_num = atoi(optarg);
            select_required = true;
            break;
        case 'h':
            PrintHelper();
            exit(EXIT_SUCCESS);
        }
    }
    CheckInputError(optind == argc);
    kmer_count_path = argv[optind++];

    if (!sample_info_path.empty())
    {
        class_num = MetaDataParser(sample_info, sample_info_path);
    }

    MethodParser(user_method_name, eval_func, fold_num, out_header, eval_method, sort_mode, class_num);
    PrintParameterInfo(kmer_count_path, sample_info_path, eval_method, fold_num, sort_mode, top_num);

    KmerEvaluate(kmer_model_vec, kmer_count_path, sample_info, fold_num, eval_func, user_method_name);

    KmerSort(kmer_model_vec, sort_mode);

    if (eval_method == METHOD_TTEST)
    {
        PValue_Adjustment_BH(kmer_model_vec);
    }

    KmerPrint(kmer_model_vec, top_num, out_header);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}