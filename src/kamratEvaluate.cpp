#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>

#include "utils/utils.hpp"
#include "utils/statistics.hpp"
#include "utils/column_info.hpp"
#include "utils/seq_coding.hpp"
#include "evaluate/parse_opt_print_info.hpp"
#include "evaluate/kmer_count_tab.hpp"

void ScanCountTable(KMerCountTab &kmer_count_tab,
                    const std::string &kmer_count_path,
                    const std::string &sample_name_path,
                    const bool stranded)
{
    std::ifstream kmer_count_file(kmer_count_path);
    ExitIf(!kmer_count_file.is_open(), "ERROR: k-mer count file " + kmer_count_path + " was not found.");
    ColumnInfo column_info;
    std::string line;
    // dealing with the header line for parsing column names //
    std::getline(kmer_count_file, line);
    column_info.MakeColumnInfo(line, sample_name_path, "NULL");
    // dealing with the count lines //
    while (std::getline(kmer_count_file, line))
    {
        std::istringstream conv(line);
        std::string term, seq;
        std::vector<float> sample_counts;
        for (size_t i(0); conv >> term; ++i)
        {
            char col_nat = column_info.GetColumnNature(i);
            if (col_nat == 'f')
            {
                seq = term;
            }
            else if (col_nat == 's')
            {
                sample_counts.push_back(std::stof(term));
            }
        }
        ExitIf(!kmer_count_tab.AddKMerCountInMem(Seq2Int(seq, seq.size(), stranded), sample_counts), "ERROR: duplicated k-mer " + seq);
    }
    kmer_count_file.close();
}

void EvalauteContigList(std::ostream &out_s,
                        const std::string &contig_list_path,
                        const KMerCountTab &kmer_count_tab,
                        const std::string &eval_method,
                        const unsigned int k_length,
                        const unsigned int min_overlap,
                        const bool stranded)
{
    std::ifstream contig_list_file(contig_list_path);
    ExitIf(!contig_list_file.is_open(), "ERROR: contig list file " + contig_list_path + "was not found.");

    std::string seq_name, seq, line;
    size_t nline(0);
    while (std::getline(contig_list_file, line))
    {
        if (nline == 0)
        {
            seq_name = line;
        }
        else if (line[0] == '>' || contig_list_file.eof())
        {
            float score = (eval_method == "mac") ? INT_MIN : INT_MAX;
            size_t left_pos = 0, right_pos = 0, nb_kmer = seq.size() - k_length + 1;
            while (right_pos < nb_kmer)
            {
                std::string left_kmer = seq.substr(left_pos, k_length), right_kmer;
                uint64_t left_code = Seq2Int(left_kmer, k_length, stranded), right_code;
                ExitIf(!kmer_count_tab.CheckKMerCodeExist(left_code), "ERROR: left k-mer " + left_kmer + " for contig " + seq + " not found in k-mer count table.");
                do
                {
                    ++right_pos;
                    right_kmer = seq.substr(right_pos, k_length);
                    right_code = Seq2Int(right_kmer, k_length, stranded);
                } while (right_pos < nb_kmer && right_pos - left_pos <= min_overlap && !kmer_count_tab.CheckKMerCodeExist(right_code));
                ExitIf(right_pos >= nb_kmer || right_pos - left_pos > min_overlap, "ERROR: no right k-mer found for " + left_kmer + " in contig " + seq + ".");
                std::vector<float> left_counts = kmer_count_tab.GetCountByKmer(left_code),
                                   right_counts = kmer_count_tab.GetCountByKmer(right_code);
                if (eval_method == "mac")
                {
                    score = GetMaxInPair<decltype(score)>(CalcMeanAbsoluteContrast(left_counts, right_counts), score);
                }
                else if (eval_method == "pearson")
                {
                    score = GetMinInPair<decltype(score)>(CalcPearsonCorrelation(left_counts, right_counts), score);
                }
                else if (eval_method == "spearman")
                {
                    score = GetMinInPair<decltype(score)>(CalcSpearmanCorrelation(left_counts, right_counts), score);
                }
                else
                {
                    ExitIf(true, "ERROR: unknown evaluate method " + eval_method);
                }
            }
            out_s << seq_name << "\t" << score << std::endl;
            seq_name = line;
            seq.clear();
        }
        else
        {
            seq += line;
        }
    }
    contig_list_file.close();
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string contig_list_path, sample_name_path, eval_method, kmer_count_path;
    unsigned int k_length(31), min_overlap(15);
    bool stranded(true);

    ParseOptions(argc, argv, contig_list_path, eval_method, stranded, k_length, min_overlap, sample_name_path, kmer_count_path);
    PrintRunInfo(contig_list_path, k_length, min_overlap, stranded, sample_name_path, eval_method, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    KMerCountTab kmer_count_tab;
    ScanCountTable(kmer_count_tab, kmer_count_path, sample_name_path, stranded);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    EvalauteContigList(std::cout, contig_list_path, kmer_count_tab, eval_method, k_length, min_overlap, stranded);

    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
