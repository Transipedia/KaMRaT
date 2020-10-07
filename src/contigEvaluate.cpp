/* ======================================================================================================================= *\
    kamratEvaluate
    Author: Haoliang Xue
  ----------------------------------------------
    For a given list of sequence, kamratEvaluate evaluates their integraty, i.e., how much are we sure that these contigs
    do come from a correct assembly/extension
    Three evaluation methods are applied: 
        1. pearson distance = 0.5 * (1 - pearson correlation)
        2. spearman distance = 0.5 * (1 - spearman correlation)
        3. distance of mean absolute contrast
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

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
    size_t pos = kmer_count_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && kmer_count_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line;
    // dealing with the header line for parsing column names
    std::getline(kmer_count_instream, line);
    kmer_count_tab.MakeSmpCond(colname_list_path);
    kmer_count_tab.MakeColumnInfo(line, "NULLFORNOSCORE");
    // dealing with the count lines
    const std::string idx_path = kmer_count_tab.GetIndexPath();
    std::ofstream idx_file;
    if (!idx_path.empty())
    {
        idx_file.open(idx_path);
        if (!idx_file.is_open()) // to assure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    for (size_t i(0); std::getline(kmer_count_instream, line); ++i)
    {
        std::istringstream conv(line);
        std::string term, seq;
        conv >> seq; // first column is supposed to be the feature column
        if (!code2serial.insert({Seq2Int(seq, seq.size(), stranded), i}).second)
        {
            throw std::domain_error("duplicated k-mer " + seq);
        }
        float _;
        kmer_count_tab.AddRowAsFields(_, line, idx_file); // don't need row score
    }
    if (idx_file.is_open())
    {
        idx_file.close();
    }
    kmer_count_file.close();
}

void EstablishSeqListFromMultilineFasta(name2seq_t &name2seq, const std::string &contig_fasta_path)
{
    std::ifstream contig_list_file(contig_fasta_path);
    if (!contig_list_file.is_open())
    {
        throw std::domain_error("contig fasta file " + contig_fasta_path + "was not found");
    }
    std::string seq_name, seq, line;
    size_t nline(0);
    while (true)
    {
        std::getline(contig_list_file, line);
        if (nline == 0) // reading the first line
        {
            seq_name = line.substr(1);
        }
        else if (line[0] == '>') // reading in middle of the file
        {
            name2seq.insert({seq_name, SeqElem(seq, nline, 0)});
            seq_name = line.substr(1);
            seq.clear();
        }
        else if (contig_list_file.eof()) // reading the last line
        {
            name2seq.insert({seq_name, SeqElem(seq, nline, 0)});
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

const void EvaluatePrintSeqElemFarthest(const std::string &name,
                                        const std::string &seq,
                                        const std::string &eval_method,
                                        const bool stranded,
                                        const unsigned int k_len,
                                        const CountTab &kmer_count_tab,
                                        const code2serial_t &code2serial,
                                        const std::string &idx_path)
{
    std::ifstream idx_file(idx_path);

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
        std::cout << name << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
    }
    else
    {
        auto iter1 = code2serial.find(Seq2Int(kmer1, k_len, stranded)),
             iter2 = code2serial.find(Seq2Int(kmer2, k_len, stranded));
        std::cout << name << "\t" << kmer_count_tab.CalcCountDistance(iter1->second, iter2->second, eval_method, idx_file) << "\t"
                  << kmer1 << "\t" << kmer2 << std::endl;
    }
    idx_file.close();
}

const void EvaluatePrintSeqElemWorstAdj(const std::string &name,
                                        const std::string &seq,
                                        const std::string &eval_method,
                                        const bool stranded,
                                        const unsigned int k_len,
                                        const CountTab &kmer_count_tab,
                                        const code2serial_t &code2serial,
                                        const std::string &idx_path)
{
    std::ifstream idx_file(idx_path);

    size_t seq_size = seq.size(), start_pos1 = 0, start_pos2;
    std::string kmer1, kmer2;
    float dist = MIN_DISTANCE - 1;
    while (start_pos1 + k_len < seq_size)
    {
        float dist_x = dist;
        std::string kmer1_x = seq.substr(start_pos1, k_len);
        while (start_pos1 + k_len < seq_size && code2serial.find(Seq2Int(kmer1_x, k_len, stranded)) == code2serial.cend())
        {
            ++start_pos1;
            if (start_pos1 + k_len < seq_size)
            {
                kmer1_x.erase(0, 1);
                kmer1_x.push_back(seq[start_pos1 + k_len - 1]);
            }
        }
        start_pos2 = start_pos1 + 1;
        std::string kmer2_x = seq.substr(start_pos2, k_len);
        while (start_pos2 + k_len <= seq_size && code2serial.find(Seq2Int(kmer2_x, k_len, stranded)) == code2serial.cend())
        {
            ++start_pos2;
            if (start_pos2 + k_len <= seq_size)
            {
                kmer2_x.erase(0, 1);
                kmer2_x.push_back(seq[start_pos2 + k_len - 1]);
            }
        }
        if (start_pos2 + k_len <= seq_size)
        {
            auto iter1 = code2serial.find(Seq2Int(kmer1_x, k_len, stranded)),
                 iter2 = code2serial.find(Seq2Int(kmer2_x, k_len, stranded));
            dist_x = kmer_count_tab.CalcCountDistance(iter1->second, iter2->second, eval_method, idx_file);
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
        std::cout << name << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
    }
    else
    {
        std::cout << name << "\t" << dist << "\t" << kmer1 << "\t" << kmer2 << std::endl;
    }

    idx_file.close();
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string contig_fasta_path, eval_method, eval_mode, colname_list_path, idx_path, kmer_count_path;
    unsigned int k_len(31);
    bool stranded(true);

    ParseOptions(argc, argv, contig_fasta_path, eval_method, eval_mode, idx_path, stranded, k_len, colname_list_path, kmer_count_path);
    PrintRunInfo(contig_fasta_path, eval_method, eval_mode, idx_path, stranded, k_len, colname_list_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    CountTab kmer_count_tab(idx_path);
    code2serial_t code2serial;
    ScanCountTable(kmer_count_tab, code2serial, kmer_count_path, colname_list_path, stranded);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    name2seq_t name2seq;
    EstablishSeqListFromMultilineFasta(name2seq, contig_fasta_path);
    std::cerr << "Number of sequence for evaluation: " << name2seq.size() << std::endl;

    std::cout << "name\teval_dist\tkmer1\tkmer2" << std::endl;
    for (const auto &ns_pair : name2seq)
    {
        auto name = ns_pair.first, seq = ns_pair.second.GetSeq();
        if (seq.size() == k_len && code2serial.find(Seq2Int(seq, k_len, stranded)) != code2serial.cend()) // if seq is a k-mer and is in k-mer count table
        {
            std::cout << name << "\t" << MIN_DISTANCE << "\t" << seq << "\t" << seq << std::endl;
        }
        else if (seq.size() <= k_len) // if seq is a k-mer but not in k-mer count table or if it's shorter than a k-mer
        {
            std::cout << name << "\t" << MAX_DISTANCE << "\tNONE\tNONE" << std::endl;
        }
        else if (eval_mode == "farthest") // if seq has multiple k-mers, and evaluation mode is "farthest"
        {
            EvaluatePrintSeqElemFarthest(name, seq, eval_method, stranded, k_len, kmer_count_tab, code2serial, idx_path);
        }
        else if (eval_mode == "worstAdj") // if seq has multiple k-mers, and evaluation mode is "worstAdj"
        {
            EvaluatePrintSeqElemWorstAdj(name, seq, eval_method, stranded, k_len, kmer_count_tab, code2serial, idx_path);
        }
    }

    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
