/* ======================================================================================================================= *\
    kamratEvaluate
    Author: Haoliang Xue
  ----------------------------------------------
    For a given list of sequence, kamratEvaluate evaluates their integraty, i.e., how much are we sure that these contigs
    do come from a correct assembly/extension
    Three evaluation methods are applied: 
        1. pearson correlation, the higher the better
        2. spearman correlation, the higher the better
        3. mean absolute contrast, the lower the better
    Two evaluation modes are applied: 
        1. by comparing for each contig the two k-mers at both ends, which are the two most distant k-mers
        2. by comparing for each contig the adjacent k-mer pairs one by one, and take the worst comparison result as the
        evaluation of the contig
    Input: 
        1. Contig list for evaluation
        2. k-mer count table
    Output: Count evluation table of four columns
        1. Contig
        2. Value for evaluation of integraty
        3. Left k-mer in related comparison
        4. Right k-mer in related comparison
\* ======================================================================================================================= */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <ctime>

#include "utils/utils.hpp"
#include "utils/statistics.hpp"
#include "utils/column_info.hpp"
#include "utils/seq_coding.hpp"
#include "data_struct/kmer_count_tab_tmplt.hpp"
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/evaluate.hpp"

const float MIN_DISTANCE = 0, MAX_DISTANCE = 1 - MIN_DISTANCE;

template <typename countT>
void ScanCountTable(KMerCountTab<countT> &kmer_count_tab,
                    const std::string &kmer_count_path,
                    const std::string &colname_list_path,
                    const bool stranded)
{
    std::ifstream kmer_count_file(kmer_count_path);
    ExitIf(!kmer_count_file.is_open(), "ERROR: k-mer count file " + kmer_count_path + " was not found");
    ColumnInfo column_info;
    std::string line;
    // dealing with the header line for parsing column names
    std::getline(kmer_count_file, line);
    column_info.MakeColumnInfo(line, colname_list_path, "NULL");
    // dealing with the count lines
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
        ExitIf(!kmer_count_tab.AddKMerCountInMem(Seq2Int(seq, seq.size(), stranded),
                                                 sample_counts),
               "ERROR: duplicated k-mer " + seq);
    }
    kmer_count_file.close();
}

void EstablishSeqListFromMultilineFasta(seqVect_t &seq_vect,
                                        const std::string &contig_list_path)
{
    std::ifstream contig_list_file(contig_list_path);
    ExitIf(!contig_list_file.is_open(), "ERROR: contig list file " + contig_list_path + "was not found.");

    std::string seq_name, seq, line;
    size_t nline(0);
    while (std::getline(contig_list_file, line))
    {
        if (nline == 0) // reading the first line
        {
            seq_name = line;
        }
        else if (line[0] == '>') // reading in middle of the file
        {
            seq_vect.emplace_back(seq_name, seq);
            seq_name = line;
            seq.clear();
        }
        else if (contig_list_file.eof()) // reading the last line
        {
            seq_vect.emplace_back(seq_name, seq);
            std::cerr << "The last line visited: " << seq << std::endl;
            seq.clear();
        }
        else
        {
            seq += line;
        }
    }
    contig_list_file.close();
}

template <typename countT>
const void EvaluatePrintSeqElemFarthest(const std::string &seq,
                                        const std::string &eval_method,
                                        const bool stranded,
                                        const unsigned int k_len,
                                        const KMerCountTab<countT> &kmer_count_tab)
{
    size_t start_pos1 = 0, start_pos2 = seq.size() - k_len;
    std::string kmer1 = seq.substr(start_pos1, k_len), kmer2 = seq.substr(start_pos2, k_len);
    std::vector<countT> kmer1_count, kmer2_count;
    while (start_pos1 < start_pos2 && !kmer_count_tab.GetCountInMem(kmer1_count, Seq2Int(kmer1, k_len, stranded)))
    {
        ++start_pos1;
	kmer1 = seq.substr(start_pos1, k_len);
    }
    while (start_pos1 < start_pos2 && !kmer_count_tab.GetCountInMem(kmer2_count, Seq2Int(kmer2, k_len, stranded)))
    {
        --start_pos2;
	kmer2 = seq.substr(start_pos2, k_len);
    }
    if (start_pos1 >= start_pos2)
    {
        std::cout << seq << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
    }
    else
    {
        std::cout << seq << "\t" << CalcDistance(kmer1_count, kmer2_count, eval_method) << "\t" << kmer1 << "\t" << kmer2;
    }
}

template <typename countT>
const void EvaluatePrintSeqElemWorstShortest(const std::string &seq,
                                             const std::string &eval_method,
                                             const bool stranded,
                                             const unsigned int k_len,
                                             const KMerCountTab<countT> &kmer_count_tab)
{
    std::string kmer1, kmer2;
    std::vector<countT> kmer1_count, kmer2_count;
    float dist = MIN_DISTANCE;
    size_t start_pos1 = 0, start_pos2 = 0;
    while (start_pos2 < seq.size())
    {
        float dist_x;
        std::string kmer1_x = seq.substr(start_pos1, k_len), kmer2_x = seq.substr(start_pos2, k_len);
        while (start_pos1 < seq.size() && !kmer_count_tab.GetCountInMem(kmer1_count, Seq2Int(kmer1_x, k_len, stranded)))
        {
            ++start_pos1;
	    kmer1_x = seq.substr(start_pos1, k_len); 
        }
        start_pos2 = start_pos1 + 1;
        while (start_pos2 < seq.size() && !kmer_count_tab.GetCountInMem(kmer2_count, Seq2Int(kmer2_x, k_len, stranded)))
        {
            ++start_pos2;
	    kmer2_x = seq.substr(start_pos2, k_len);
        }
        start_pos1 = start_pos2;
        dist_x = CalcDistance(kmer1_count, kmer2_count, eval_method);
        if (dist_x > dist)
        {
            dist = dist_x;
            kmer1 = kmer1_x;
            kmer2 = kmer2_x;
        }
    }
    std::cout << seq << "\t" << dist << "\t" << kmer1 << "\t" << kmer2 << std::endl;
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string contig_list_path, eval_method, eval_mode, colname_list_path, kmer_count_path;
    unsigned int k_len(31);
    bool stranded(true);

    ParseOptions(argc, argv, contig_list_path, eval_method, eval_mode, stranded, k_len, colname_list_path, kmer_count_path);
    PrintRunInfo(contig_list_path, eval_method, eval_mode, stranded, k_len, colname_list_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    KMerCountTab<float> kmer_count_tab;
    ScanCountTable(kmer_count_tab, kmer_count_path, colname_list_path, stranded);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    seqVect_t seq_vect;
    EstablishSeqListFromMultilineFasta(seq_vect, contig_list_path);

    std::cout << "sequence\teval_dist\tkmer1\tkmer2" << std::endl;
    for (const auto &seq_elem : seq_vect)
    {
        auto seq = seq_elem.GetSeq();
        if (seq.size() == k_len && kmer_count_tab.IsKMerExist(Seq2Int(seq, k_len, stranded))) // if seq has only one k-mer and is in k-mer count table
        {
            std::cout << seq << "\t" << MIN_DISTANCE << "\t" << seq << "\t" << seq << std::endl;
        }
        else if (seq.size() == k_len) // if seq has only one k-mer but is not in k-mer count table
        {
            std::cout << seq << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
        }
        else if (eval_mode == "farthest") // if seq has multiple k-mers, and evaluation mode is "farthest"
        {
            EvaluatePrintSeqElemFarthest(seq, eval_method, stranded, k_len, kmer_count_tab);
        }
        else
        {
            EvaluatePrintSeqElemWorstShortest(seq, eval_method, stranded, k_len, kmer_count_tab);
        }
    }

    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
