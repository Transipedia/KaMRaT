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
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>

#include "utils/statistics.hpp"
#include "utils/seq_coding.hpp"
#include "data_struct/count_tab.hpp"
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/evaluate.hpp"

const float MIN_DISTANCE = 0, MAX_DISTANCE = 1 - MIN_DISTANCE;

void ScanCountTable(CountTab &kmer_count_tab,
                        code2serial_t &code2serial,
                        const std::string &kmer_count_path,
                        const std::string &colname_list_path,
                        const bool stranded)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        throw std::domain_error("k-mer count file " + kmer_count_path + " was not found");
    }
    std::string line;
    // dealing with the header line for parsing column names
    std::getline(kmer_count_file, line);
    kmer_count_tab.MakeColumnInfo(line, colname_list_path, "NULLFORNOSCORE");
    // dealing with the count lines
    while (std::getline(kmer_count_file, line))
    {
        std::istringstream conv(line);
        std::string term, seq;
        conv >> seq; // first column is supposed to be the feature column
        if (!code2serial.insert({Seq2Int(seq, seq.size(), stranded), static_cast<size_t>(kmer_count_tab.AddCountInMem(line) + 0.5)}).second)
        {
            throw std::domain_error("duplicated k-mer " + seq);
        }
    }
    kmer_count_file.close();
}

void EstablishSeqListFromMultilineFasta(tag2seq_t &tag_seq_dict,
                                        const std::string &contig_fasta_path)
{
    std::ifstream contig_list_file(contig_fasta_path);
    if (!contig_list_file.is_open())
    {
        throw std::domain_error("contig fasta file " + contig_fasta_path + "was not found");
    }
    std::string seq_tag, seq, line;
    size_t nline(0);
    while (true)
    {
        std::getline(contig_list_file, line);
        if (nline == 0) // reading the first line
        {
            seq_tag = line.substr(1);
        }
        else if (line[0] == '>') // reading in middle of the file
        {
            tag_seq_dict.insert({seq_tag, SeqElem(seq)});
            seq_tag = line.substr(1);
            seq.clear();
        }
        else if (contig_list_file.eof()) // reading the last line
        {
            tag_seq_dict.insert({seq_tag, SeqElem(seq)});
            seq.clear();
            break;
        }
        else
        {
            seq += line;
        }
        ++nline;
    }
    contig_list_file.close();
}

const void EvaluatePrintSeqElemFarthest(const std::string &tag,
                                        const std::string &seq,
                                        const std::string &eval_method,
                                        const bool stranded,
                                        const unsigned int k_len,
                                        const CountTab &kmer_count_tab,
                                        const code2serial_t &code2serial)
{
    size_t start_pos1 = 0, start_pos2 = seq.size() - k_len;
    std::string kmer1 = seq.substr(start_pos1, k_len), kmer2 = seq.substr(start_pos2, k_len);
    while (start_pos1 < start_pos2 && code2serial.find(Seq2Int(kmer1, k_len, stranded)) == code2serial.cend())
    {
        ++start_pos1;
        kmer1 = seq.substr(start_pos1, k_len);
    }
    while (start_pos1 < start_pos2 && code2serial.find(Seq2Int(kmer2, k_len, stranded)) == code2serial.cend())
    {
        --start_pos2;
        kmer2 = seq.substr(start_pos2, k_len);
    }
    if (start_pos1 >= start_pos2)
    {
        std::cout << tag << "\t" << seq << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
    }
    else
    {
        auto iter1 = code2serial.find(Seq2Int(kmer1, k_len, stranded)),
             iter2 = code2serial.find(Seq2Int(kmer2, k_len, stranded));
        std::vector<float> kmer1_count, kmer2_count;
        kmer_count_tab.GetCountInMem(kmer1_count, iter1->second);
        kmer_count_tab.GetCountInMem(kmer2_count, iter2->second);
        std::cout << tag << "\t" << seq << "\t" << CalcDistance(kmer1_count, kmer2_count, eval_method) << "\t" << kmer1 << "\t" << kmer2 << std::endl;
    }
}

const void EvaluatePrintSeqElemWorstAdj(const std::string &tag,
                                        const std::string &seq,
                                        const std::string &eval_method,
                                        const bool stranded,
                                        const unsigned int k_len,
                                        const CountTab &kmer_count_tab,
                                        const code2serial_t &code2serial)
{
    std::string kmer1, kmer2;
    float dist = MIN_DISTANCE - 1;
    size_t start_pos1 = 0, start_pos2;
    while (start_pos1 + k_len < seq.size())
    {
        float dist_x = dist;
        std::string kmer1_x = seq.substr(start_pos1, k_len);
        while (start_pos1 + k_len <= seq.size() && code2serial.find(Seq2Int(kmer1_x, k_len, stranded)) == code2serial.cend())
        {
            ++start_pos1;
            kmer1_x.erase(0, 1);
            kmer1_x.push_back(seq.at(start_pos1 + k_len - 1));
        }
        start_pos2 = start_pos1 + 1;
        std::string kmer2_x = seq.substr(start_pos2, k_len);
        while (start_pos2 + k_len <= seq.size() && code2serial.find(Seq2Int(kmer2_x, k_len, stranded)) == code2serial.cend())
        {
            ++start_pos2;
            kmer2_x.erase(0, 1);
            kmer2_x.push_back(seq.at(start_pos2 + k_len - 1));
        }
        if (start_pos2 + k_len <= seq.size())
        {
            auto iter1 = code2serial.find(Seq2Int(kmer1_x, k_len, stranded)),
                 iter2 = code2serial.find(Seq2Int(kmer2_x, k_len, stranded));
            std::vector<float> kmer1_count, kmer2_count;
            kmer_count_tab.GetCountInMem(kmer1_count, iter1->second);
            kmer_count_tab.GetCountInMem(kmer2_count, iter2->second);
            dist_x = CalcDistance(kmer1_count, kmer2_count, eval_method);
            if (dist_x > dist)
            {
                dist = dist_x;
                kmer1 = kmer1_x;
                kmer2 = kmer2_x;
            }
        }
        start_pos1 = start_pos2;
    }
    if (kmer1.empty() || kmer2.empty())
    {
        std::cout << tag << "\t" << seq << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
    }
    else
    {
        std::cout << tag << "\t" << seq << "\t" << dist << "\t" << kmer1 << "\t" << kmer2 << std::endl;
    }
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string contig_fasta_path, eval_method, eval_mode, colname_list_path, kmer_count_path;
    unsigned int k_len(31);
    bool stranded(true);

    if (!ParseOptions(argc, argv, contig_fasta_path, eval_method, eval_mode, stranded, k_len, colname_list_path, kmer_count_path))
    {
        throw std::domain_error("option parsing interrupted, please check kamratEvaluate -h");
    }
    PrintRunInfo(contig_fasta_path, eval_method, eval_mode, stranded, k_len, colname_list_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    CountTab kmer_count_tab("inMem");
    code2serial_t code2serial;
    ScanCountTable(kmer_count_tab, code2serial, kmer_count_path, colname_list_path, stranded);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    tag2seq_t tag_seq_dict;
    EstablishSeqListFromMultilineFasta(tag_seq_dict, contig_fasta_path);
    std::cerr << "Number of sequence for evaluation: " << tag_seq_dict.size() << std::endl;

    std::cout << "tag\tseq\teval_dist\tkmer1\tkmer2" << std::endl;
    for (const auto &ns_pair : tag_seq_dict)
    {
        auto tag = ns_pair.first, seq = ns_pair.second.GetSeq();
        if (seq.size() == k_len && code2serial.find(Seq2Int(seq, k_len, stranded)) != code2serial.cend()) // if seq has only one k-mer and is in k-mer count table
        {
            std::cout << tag << "\t" << seq << "\t" << MIN_DISTANCE << "\t" << seq << "\t" << seq << std::endl;
        }
        else if (seq.size() == k_len) // if seq has only one k-mer but is not in k-mer count table
        {
            std::cout << tag << "\t" << seq << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
        }
        else if (eval_mode == "farthest") // if seq has multiple k-mers, and evaluation mode is "farthest"
        {
            EvaluatePrintSeqElemFarthest(tag, seq, eval_method, stranded, k_len, kmer_count_tab, code2serial);
        }
        else if (eval_mode == "worstAdj")
        {
            EvaluatePrintSeqElemWorstAdj(tag, seq, eval_method, stranded, k_len, kmer_count_tab, code2serial);
        }
    }

    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
