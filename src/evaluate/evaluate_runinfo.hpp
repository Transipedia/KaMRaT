#ifndef KAMRAT_RUNINFOPARSER_EVALUATE_HPP
#define KAMRAT_RUNINFOPARSER_EVALUATE_HPP

#include <iostream>
#include <string>
#include <set>

const std::set<std::string> kEvaluateMethodUniv{"pearson", "spearman", "mac", "mean", "median"};

using fastaVect_t = std::vector<std::pair<std::string, std::string>>;

inline void PrintContigEvaluateHelper()
{
    std::cerr << "========= kamratEvaluate helper =========" << std::endl;
    std::cerr << "[Usage]    contigEvaluate -klen INT -fasta STR -eval-method STR -idx-path STR [-options] KMER_COUNT_TAB" << std::endl
              << std::endl;
    std::cerr << "[Option]    -h, -help           Print the helper" << std::endl;
    std::cerr << "            -klen INT           k-mer length, MANDATORY" << std::endl;
    std::cerr << "            -fasta STR          Contig fasta file path, MANDATORY" << std::endl;
    std::cerr << "            -eval-method STR    Evaluate method, MANDATORY, method can be:" << std::endl
              << "                                    pearson     evaluate contig's largest Pearson distance between findable adjacent k-mers" << std::endl
              << "                                    spearman    evaluate contig's largest Spearman distance between findable adjacent k-mers" << std::endl
              << "                                    mac         evaluate contig's largest MAC distance between findable adjacent k-mers" << std::endl
              << "                                    mean        evaluate contig's mean count among all composite k-mers for each sample" << std::endl
              << "                                    median      evaluate contig's median count among all composite k-mers for each sample" << std::endl;
    std::cerr << "            -idx-path STR       Count index file path, MANDATORY" << std::endl;
    std::cerr << "            -unstrand           If the k-mers are generated from non-stranded RNA-seq data" << std::endl;
    std::cerr << "            -smp-info STR       Sample-info path, either list or table with sample names as the first column" << std::endl
              << "                                    if absent, all columns except the first one in k-mer count table are taken as samples" << std::endl;
    std::cerr << "            -no-norm            Evaluate WITHOUT count nomralization" << std::endl;
    std::cerr << "            -contig-name        Output contig names instead of contig sequence" << std::endl;
    std::cerr << "            -max-shift INT      Maximum allowed shift between k-mers [default: inf]" << std::endl;
    std::cerr << "            -out-path STR       Output table path [default: output to screen]" << std::endl
              << std::endl;
    std::cerr << "[Output]    if evaluation method is: pearson, spearman, or mac, the output will be 4 columns for each contig:" << std::endl
              << "                contig          contig name or sequence" << std::endl
              << "                dist            evaluated distance by the given method" << std::endl
              << "                kmer1           left k-mer in the pair that gives the largest distance" << std::endl
              << "                kmer2           right k-mer in the pair that gives the largest distance" << std::endl;
    std::cerr << "            if evaluation method is: mean or median, the output will be a matrix:" << std::endl
              << "                first column is the contig name or sequence" << std::endl
              << "                next columns indicate sample counts" << std::endl;
}

inline void PrintRunInfo(const size_t k_len,
                         const std::string &contig_fasta_path,
                         const std::string &eval_method,
                         const std::string &idx_path,
                         const bool stranded,
                         const std::string &colname_list_path,
                         const bool no_norm,
                         const bool contig_name,
                         const std::string &out_path,
                         const std::string &kmer_count_path)
{
    std::cerr << std::endl;
    std::cerr << "k-mer length:            " << k_len << std::endl;
    std::cerr << "Contig list path:        " << contig_fasta_path << std::endl;
    std::cerr << "Evaluation method:       " << eval_method << std::endl;
    std::cerr << "Index directory path:    " << idx_path << std::endl;
    std::cerr << "Stranded mode:           " << (stranded ? "On" : "Off") << std::endl;
    if (!colname_list_path.empty())
    {
        std::cerr << "Colname list path:       " << colname_list_path << std::endl;
    }
    std::cerr << "Count normalization:     " << (no_norm ? "no" : "yes") << std::endl;
    std::cerr << "Contig output:           " << (contig_name ? "by name" : "by sequence") << std::endl;
    if (!out_path.empty())
    {
        std::cerr << "Output path:             " << out_path << std::endl;
    }
    else
    {
        std::cerr << "Output to screen" << std::endl;
    }
    std::cerr << "k-mer count path:        " << kmer_count_path << std::endl;
}

inline void ParseOptions(int argc,
                         char *argv[],
                         size_t &k_len,
                         std::string &contig_fasta_path,
                         std::string &eval_method,
                         std::string &idx_path,
                         bool &stranded,
                         std::string &colname_list_path,
                         bool &no_norm,
                         bool &contig_name,
                         std::string &out_path,
                         std::string &kmer_count_path)
{
    int i_opt = 1;
    while (i_opt < argc && argv[i_opt][0] == '-')
    {
        std::string arg(argv[i_opt]);
        if (arg == "-help" || arg == "-h")
        {
            PrintContigEvaluateHelper();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-klen" && i_opt + 1 < argc)
        {
            k_len = atoi(argv[++i_opt]);
        }
        else if (arg == "-fasta" && i_opt + 1 < argc)
        {
            contig_fasta_path = argv[++i_opt];
        }
        else if (arg == "-eval-method" && i_opt + 1 < argc)
        {
            eval_method = argv[++i_opt];
        }
        else if (arg == "-idx-path" && i_opt + 1 < argc)
        {
            idx_path = argv[++i_opt];
        }
        else if (arg == "-unstrand")
        {
            stranded = false;
        }
        else if (arg == "-smp-info" && i_opt + 1 < argc)
        {
            colname_list_path = argv[++i_opt];
        }
        else if (arg == "-no-norm")
        {
            no_norm = true;
        }
        else if (arg == "-contig-name") 
        {
            contig_name = true;
        }
        else if (arg == "-out-path" && i_opt + 1 < argc)
        {
            out_path = argv[++i_opt];
        }
        else
        {
            PrintContigEvaluateHelper();
            throw std::invalid_argument("unknown option " + std::string(argv[i_opt]));
        }
        ++i_opt;
    }
    if (i_opt == argc)
    {
        PrintContigEvaluateHelper();
        throw std::invalid_argument("k-mer count table path is mandatory");
    }
    kmer_count_path = argv[i_opt++];

    if (k_len == 0)
    {
        PrintContigEvaluateHelper();
        throw std::invalid_argument("k-mer length is missed or failed to be parsed");
    }
    if (contig_fasta_path.empty())
    {
        PrintContigEvaluateHelper();
        throw std::invalid_argument("contig path is missed or failed to be parsed");
    }
    if (kEvaluateMethodUniv.find(eval_method) == kEvaluateMethodUniv.cend())
    {
        PrintContigEvaluateHelper();
        throw std::invalid_argument("evaluate method is missed or failed to be parsed: " + eval_method);
    }
    if (idx_path.empty())
    {
        PrintContigEvaluateHelper();
        throw std::invalid_argument("Index path is missed or failed to be parsed");
    }
}

#endif //KAMRAT_RUNINFOPARSER_EVALUATE_HPP
