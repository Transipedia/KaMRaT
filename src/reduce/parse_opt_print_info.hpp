#ifndef KAMRAT_REDUCE_PARSEOPTPRINTINFO_HPP
#define KAMRAT_REDUCE_PARSEOPTPRINTINFO_HPP

#include <iostream>
#include <string>

#include "utils.hpp"
#include "eval_methods.hpp"

inline void PrintHelper()
{
    std::cerr << "=====> kamratReduce Helper <=====" << std::endl;
    std::cerr << "[Usage]    kamratReduce [-d samp_info_path -m eval_method:fold_num -s sort_mode -N top_num -T transf_mode -C count_offset] kmer_count_path" << std::endl
              << std::endl;
    std::cerr << "[Option]        -d STRING    Path to sample-condition or sample file, without header line" << std::endl
              << "                             if absent, all except the first column in k-mer count table will be regarded as samples" << std::endl;
    std::cerr << "                -m STRING    Evaluation method name and fold number for cross-validation (if needed), seperated by \':\'" << std::endl;
    std::cerr << "                -s STRING    Sorting mode, default value depends on evaluation method (c.f [SORT MODE])" << std::endl;
    std::cerr << "                -N INT       Number of top features to print" << std::endl;
    std::cerr << "                -T STRING    Transformation before evaluation (e.g. log)" << std::endl;
    std::cerr << "                -C DOUBLE    Count offset to be added before transformation and evaluation [0 without log, 1 with log transformation or log2FC evaluation method]" << std::endl;
    std::cerr << "                -r INT       Accepted minimum recurrence number among samples" << std::endl;
    std::cerr << "                -a INT       Accepted minimum count abundance for counting sample recurrence" << std::endl
              << std::endl;
    std::cerr << "[Evaluation]    nb           Naive Bayes classification between conditions" << std::endl;
    std::cerr << "                lr           Logistic regression (slower than Naive Bayes) between conditions" << std::endl;
    std::cerr << "                sd           Standard deviation" << std::endl;
    std::cerr << "                rsd          Relative standard deviation" << std::endl;
    std::cerr << "                mc           Contrast of mean between conditions" << std::endl;
    std::cerr << "                rsdc         Contrast of relative standard deviation between conditions" << std::endl;
    std::cerr << "                ttest        T-test between conditions" << std::endl;
    std::cerr << "                es           Effect size between conditions" << std::endl;
    std::cerr << "                lfc:mean     Log2 fold change by group mean, 'mean' can be omitted as default value" << std::endl;
    std::cerr << "                lfc:median   Log2 fold change by group median" << std::endl;
    std::cerr << "                user:name    User-defined method, where name indicates a column in the k-mer count table" << std::endl
              << std::endl;
    std::cerr << "[SORT MODE]     dec          Sorting by decreasing order                              (as default for nb, lr, sd, rsd, user:name)" << std::endl;
    std::cerr << "                dec:abs      Sorting by decreasing order but on the absolute value    (as default for mdc, rsdc, es, lfc:mean, lfc:median)" << std::endl;
    std::cerr << "                inc          Sorting by increasing order                              (as default for ttest)" << std::endl;
    std::cerr << "                inc:abs      Sorting by increasing order but on the absolute value" << std::endl
              << std::endl;
}

inline void PrintRunInfo(const std::string &kmer_count_path,
                         const std::string &sample_info_path,
                         const std::string &eval_method,
                         const size_t fold_num,
                         const std::string &sort_mode,
                         const size_t top_num,
                         const std::string &transf_mode,
                         const double count_offset,
                         const size_t min_rec,
                         const size_t min_rec_abd)
{
    std::cerr << "k-mer count path:                " << kmer_count_path << std::endl;
    if (!sample_info_path.empty())
    {
        std::cerr << "Sample info path:                " << sample_info_path << std::endl;
    }
    std::cerr << "Evaluation method:               " << eval_method << std::endl;
    if (eval_method == METHOD_NB || eval_method == METHOD_LR)
    {
        std::cerr << "Fold number:                     " << fold_num << std::endl;
    }
    if (!sort_mode.empty())
    {
        std::cerr << "Sorting mode:                    " << sort_mode << std::endl;
    }
    if (top_num > 0)
    {
        std::cerr << "Number of k-mers to select:      " << top_num << std::endl;
    }
    if (!transf_mode.empty())
    {
        std::cerr << "Transformation mode:             " << transf_mode << std::endl;
    }
    if (count_offset != 0)
    {
        std::cerr << "Count offset:                    " << count_offset << std::endl;
    }
    if (min_rec != 0)
    {
        std::cerr << "Minimum recurrence:              " << min_rec << std::endl;
    }
    if (min_rec_abd != 0)
    {
        std::cerr << "Minimum recurrence abundance:    " << min_rec_abd << std::endl;
    }
    std::cerr << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         std::string &sample_info_path,
                         std::string &eval_method,
                         std::string &sort_mode,
                         size_t &top_num,
                         std::string &kmer_count_path,
                         std::string &transf_mode,
                         double &count_offset,
                         size_t &min_rec,
                         size_t &min_rec_abd,
                         bool &is_offset_set)
{
    int i_opt(1);
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-h")
        {
            PrintHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-d" && i_opt + 1 < argc)
        {
            sample_info_path = argv[++i_opt];
        }
        else if (arg == "-m" && i_opt + 1 < argc)
        {
            eval_method = argv[++i_opt];
        }
        else if (arg == "-s" && i_opt + 1 < argc)
        {
            sort_mode = argv[++i_opt];
        }
        else if (arg == "-N" && i_opt + 1 < argc)
        {
            top_num = atoi(argv[++i_opt]);
        }
        else if (arg == "-T" && i_opt + 1 < argc)
        {
            transf_mode = argv[++i_opt];
        }
        else if (arg == "-C" && i_opt + 1 < argc)
        {
            count_offset = atof(argv[++i_opt]);
            is_offset_set = true;
        }
        else if (arg == "-r" && i_opt + 1 < argc)
        {
            min_rec = atoi(argv[++i_opt]);
        }
        else if (arg == "-a" && i_opt + 1 < argc)
        {
            min_rec_abd = atoi(argv[++i_opt]);
        }
        else
        {
            std::cerr << "ERROR: unknown option " << argv[i_opt] << std::endl;
            exit(EXIT_FAILURE);
        }
        ++i_opt;
    }
    if (transf_mode == "log" && !is_offset_set)
    {
        count_offset = 1;
    }
    if (i_opt == argc)
    {
        std::cerr << "ERROR: k-mer count table path is mandatory." << std::endl
                  << "       please type kamratReduce -h for help..." << std::endl;
        exit(EXIT_FAILURE);
    }
    kmer_count_path = argv[i_opt++];
}

inline void SubCommandParser(std::string &command, std::string &sub_command)
{
    size_t split_pos = command.find(":");
    if (split_pos != std::string::npos)
    {
        sub_command = command.substr(split_pos + 1);
        command = command.substr(0, split_pos);
    }
}

void ParseMethod(std::string &user_method_name,
                 evalfuncptr_t &eval_func,
                 size_t &nb_fold,
                 std::string &out_header,
                 std::string &eval_method,
                 std::string &sort_mode,
                 double &count_offset,
                 const size_t nb_condition,
                 const bool is_offset_set)
{
    std::string sub_eval_method;

    SubCommandParser(eval_method, sub_eval_method);

    out_header = "feature\t";
    if (eval_method == METHOD_NB)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: naive Bayes evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "f1.naiveBayes";
        eval_func = nb_evaluation;
        if (!sub_eval_method.empty())
        {
            nb_fold = std::stoi(sub_eval_method);
        }
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (eval_method == METHOD_LR)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: logistic regression evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "f1.LogitReg";
        eval_func = lr_evaluation;
        if (!sub_eval_method.empty())
        {
            nb_fold = std::stoi(sub_eval_method);
        }
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (eval_method == METHOD_SD)
    {
        out_header += "sd";
        eval_func = sd_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (eval_method == METHOD_RSD)
    {
        out_header += "relat_sd";
        eval_func = rsd_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC;
        }
    }
    else if (eval_method == METHOD_MC)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: mean contrast evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "mean.ctrst";
        eval_func = mc_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC_ABS;
        }
    }
    else if (eval_method == METHOD_RSDC)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: relative standard deviation evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "rsd.ctrst";
        eval_func = rsdc_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC_ABS;
        }
    }
    else if (eval_method == METHOD_TTEST)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: t-test evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "pval.t-test";
        eval_func = ttest_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_INC;
        }
    }
    else if (eval_method == METHOD_ES)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: effect size evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        out_header += "effect_size";
        eval_func = es_evaluation;
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC_ABS;
        }
    }
    else if (eval_method == METHOD_LFC)
    {
        if (nb_condition != 2)
        {
            std::cerr << "ERROR: log fold change evaluation can only be applied on bi-condition data." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (sub_eval_method.empty() || sub_eval_method == "mean")
        {
            out_header += "Log2FC:mean";
            eval_func = lfcmean_evaluation;
            eval_method += ":mean";
        }
        else if (sub_eval_method == "median")
        {
            out_header += "Log2FC:median";
            eval_func = lfcmedian_evaluation;
            eval_method += ":median";
        }
        if (sort_mode.empty())
        {
            sort_mode = SORT_DEC_ABS;
        }
        if (!is_offset_set)
        {
            count_offset = 1;
        }
        else if (count_offset == 0)
        {
            std::cerr << "ERROR: evaluating with log2FoldChange but with count_offset=0 can be problematic." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (eval_method == METHOD_USER)
    {
        if (sub_eval_method.empty())
        {
            std::cerr << "ERROR: missing user-defined method name." << std::endl
                      << "       please type kamratReduce -h for help..." << std::endl;
            exit(EXIT_FAILURE);
        }
        user_method_name = sub_eval_method;
        out_header += user_method_name;
        eval_func = NULL;
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

void ParseSampleLabel(arma::Row<size_t> &sample_labels,
                      const std::vector<int> &header_term_nat)
{
    std::vector<int> label_vect;

    for (size_t i(0); i < header_term_nat.size(); ++i)
    {
        if (header_term_nat.at(i) >= 0)
        {
            label_vect.push_back(header_term_nat.at(i));
        }
    }
    sample_labels.set_size(label_vect.size());
    sample_labels = arma::conv_to<arma::Row<size_t>>::from(label_vect);
}

#endif //KAMRAT_REDUCE_PARSEOPTPRINTINFO_HPP