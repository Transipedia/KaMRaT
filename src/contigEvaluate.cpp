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
#include "utils/vec_operation.hpp"
#include "data_struct/tab_header.hpp"
#include "data_struct/tab_elem.hpp"
#include "run_info_parser/evaluate.hpp"

const float MIN_DISTANCE = 0, MAX_DISTANCE = 1 - MIN_DISTANCE;

void ScanCountTable(TabHeader &tab_header, featuretab_t &kmer_count_tab, code2serial_t &code2serial,
                    const size_t k_len, const bool stranded,
                    const std::string &kmer_count_path, const std::string &idx_path)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        throw std::domain_error("k-mer count file " + kmer_count_path + " was not found");
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = kmer_count_path.find_last_of(".");
        if (pos != std::string::npos && kmer_count_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line, seq;
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::getline(kmer_count_instream, line);
    tab_header.MakeColumnInfo(line, "");
    //----- Dealing with Following k-mer Count Lines -----//
    std::ofstream idx_file;
    if (!idx_path.empty())
    {
        idx_file.open(idx_path);
        if (!idx_file.is_open()) // to ensure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    std::istringstream conv;
    float _;
    for (size_t iline(0); std::getline(kmer_count_instream, line); ++iline)
    {
        conv.str(line);
        kmer_count_tab.emplace_back(conv, idx_file, _, tab_header);
        seq = std::move(line.substr(0, line.find_first_of(" \t"))); // first column as feature (string)
        if (seq.size() != k_len)
        {
            throw std::domain_error("the given k-len parameter not coherent with input k-mer: " + seq);
        }
        uint64_t kmer_uniqcode = Seq2Int(seq, seq.size(), stranded);
        if (!code2serial.insert({kmer_uniqcode, iline}).second)
        {
            throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
        }
        conv.clear();
    }
    if (idx_file.is_open())
    {
        idx_file.close();
    }
    kmer_count_tab.shrink_to_fit();
    kmer_count_file.close();
}

void EstablishSeqListFromMultilineFasta(std::unordered_map<std::string, std::string> &name2seq,
                                        const std::string &contig_fasta_path)
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
            name2seq.insert({seq_name, seq});
            seq_name = line.substr(1);
            seq.clear();
        }
        else if (contig_list_file.eof()) // reading the last line
        {
            name2seq.insert({seq_name, seq});
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

const void EvaluatePrintSeqElemFarthest(const std::string &name, const std::string &seq,
                                        const size_t k_len, const bool stranded, const std::string &eval_method,
                                        const featuretab_t &kmer_count_tab, const code2serial_t &code2serial,
                                        const std::string &idx_path, const size_t nb_count)
{
    static std::vector<float> pred_count_vect, succ_count_vect;
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
        std::cout << name << "\t"
                  << CalcXDist(kmer_count_tab[iter1->second].GetCountVect(pred_count_vect, idx_file, nb_count),
                               kmer_count_tab[iter2->second].GetCountVect(succ_count_vect, idx_file, nb_count), eval_method)
                  << "\t" << kmer1 << "\t" << kmer2 << std::endl;
    }
    idx_file.close();
}

const void EvaluatePrintSeqElemWorstAdj(const std::string &name, const std::string &seq,
                                        const size_t k_len, const bool stranded, const std::string &eval_method,
                                        const featuretab_t &kmer_count_tab, const code2serial_t &code2serial,
                                        const std::string &idx_path, const size_t nb_count)
{
    static std::vector<float> pred_count_vect, succ_count_vect;
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
            dist_x = CalcXDist(kmer_count_tab[iter1->second].GetCountVect(pred_count_vect, idx_file, nb_count),
                               kmer_count_tab[iter2->second].GetCountVect(succ_count_vect, idx_file, nb_count), eval_method);
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
    std::string contig_fasta_path, eval_method, eval_mode, colname_list_path, idx_path, kmer_count_path, out_path;
    unsigned int k_len(31);
    bool stranded(true);

    ParseOptions(argc, argv, contig_fasta_path, eval_method, eval_mode, idx_path, stranded, k_len, colname_list_path, out_path, kmer_count_path);
    PrintRunInfo(contig_fasta_path, eval_method, eval_mode, idx_path, stranded, k_len, colname_list_path, out_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    TabHeader tab_header(colname_list_path);
    featuretab_t kmer_count_tab;
    code2serial_t code2serial;
    ScanCountTable(tab_header, kmer_count_tab, code2serial, k_len, stranded, kmer_count_path, idx_path);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::unordered_map<std::string, std::string> name2seq;
    EstablishSeqListFromMultilineFasta(name2seq, contig_fasta_path);
    std::cerr << "Number of sequence for evaluation: " << name2seq.size() << std::endl;

    std::cout << "name\teval_dist\tkmer1\tkmer2" << std::endl;
    for (const auto &ns_pair : name2seq)
    {
        auto name = ns_pair.first, seq = ns_pair.second;
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
            EvaluatePrintSeqElemFarthest(name, seq, k_len, stranded, eval_method, kmer_count_tab, code2serial, idx_path, tab_header.GetNbCount());
        }
        else if (eval_mode == "worstAdj") // if seq has multiple k-mers, and evaluation mode is "worstAdj"
        {
            EvaluatePrintSeqElemWorstAdj(name, seq, k_len, stranded, eval_method, kmer_count_tab, code2serial, idx_path, tab_header.GetNbCount());
        }
    }

    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
