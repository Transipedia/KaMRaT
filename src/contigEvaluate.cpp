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

#include "common/seq_coding.hpp"
#include "common/tab_header.hpp"
#include "common/vec_operation.hpp"
#include "common/kmer_elem.hpp"
#include "evaluate/evaluate_runinfo.hpp"

const float kMinDistance = 0, kMaxDistance = 1;

void ScanCountTable(TabHeader &tab_header, code2kmer_t &code2kmer, std::vector<double> &nf_vect, const size_t k_len, const bool stranded,
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

    std::string line;
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::getline(kmer_count_instream, line);
    std::istringstream conv(line);
    tab_header.MakeColumns(conv, "");
    conv.clear();
    nf_vect.resize(tab_header.GetNbCount(), 0);
    std::cerr << "\t => Number of sample parsed: " << nf_vect.size() << std::endl;

    std::ofstream idx_file(idx_path);
    if (!idx_file.is_open()) // to ensure the file is opened
    {
        throw std::domain_error("error open file: " + idx_path);
    }

    //----- Dealing with Following k-mer Count Lines -----//
    std::vector<float> count_vect;
    std::string value_str, seq;
    for (size_t iline(0); std::getline(kmer_count_instream, line); ++iline)
    {
        conv.str(line);
        float rep_value = tab_header.ParseRowStr(count_vect, value_str, conv);
        for (size_t ismp(0); ismp < count_vect.size(); ++ismp)
        {
            nf_vect[ismp] += count_vect[ismp];
        }
        seq = line.substr(0, line.find_first_of(" \t")); // first column as feature (string)
        if (seq.size() != k_len)
        {
            throw std::domain_error("the given k-len parameter not coherent with input k-mer: " + seq);
        }
        if (!code2kmer.insert({Seq2Int(seq, seq.size(), stranded), KMerElem(rep_value, count_vect, value_str, idx_file)}).second)
        {
            throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list): " + seq);
        }
        conv.clear();
    }
    double mean_sample_sum = (std::accumulate(nf_vect.cbegin(), nf_vect.cend(), 0.0) / nf_vect.size());
    for (size_t i(0); i < nf_vect.size(); ++i)
    {
        if (nf_vect[i] == 0)
        {
            nf_vect[i] = 0;
        }
        else
        {
            nf_vect[i] = mean_sample_sum / nf_vect[i];
        }
    }
    kmer_count_file.close();
    idx_file.close();
}

void EstablishSeqListFromMultilineFasta(fastaVect_t &seq_vect, const std::string &contig_fasta_path)
{
    std::ifstream contig_list_file(contig_fasta_path);
    if (!contig_list_file.is_open())
    {
        throw std::domain_error("contig fasta file " + contig_fasta_path + "was not found");
    }
    std::string seq_name, seq, line;
    bool is_first_line(true);
    while (std::getline(contig_list_file, line))
    {
        if (line[0] == '>') // if a header line (i.e. a new sequence)
        {
            if (!is_first_line)
            {
                seq_vect.emplace_back(std::make_unique<seqElem_t>(std::make_pair(seq_name, seq))); // cache the previous sequence
                seq.clear();                                                                       // prepare for the next sequence
            }
            seq_name = line.substr(1); // record the next sequence name
            is_first_line = false;
        }
        else
        {
            seq += line;
        }
    }
    seq_vect.emplace_back(std::make_unique<seqElem_t>(std::make_pair(seq_name, seq))); // cache the last sequence
    contig_list_file.close();
}

const float EvaluateSeqDist(std::string &max_kmer1, std::string &max_kmer2,
                            const std::string &contig_seq, const std::vector<double> &nf_vect, const bool no_norm,
                            const size_t k_len, const bool stranded, const std::string &eval_method,
                            const code2kmer_t &code2kmer, std::ifstream &idx_file, const size_t nb_count, const size_t max_shift)
{
    if (contig_seq.size() == k_len && code2kmer.find(Seq2Int(contig_seq, k_len, stranded)) != code2kmer.cend()) // if seq is a k-mer in k-mer count table
    {
        max_kmer1 = max_kmer2 = contig_seq;
        return kMinDistance;
    }
    else if (contig_seq.size() <= k_len) // if seq is a k-mer but not in k-mer count table or if it's shorter than a k-mer
    {
        max_kmer1 = max_kmer2 = "NONE";
        return kMaxDistance;
    }

    static std::vector<float> left_count_vect, right_count_vect;
    size_t seq_size = contig_seq.size(), max_kmer_code1, max_kmer_code2;

    // Start condition: left k-mer as the first findable one in the k-mer count table //
    size_t start_pos1 = 0;
    uint64_t kmer_code1 = Seq2Int(contig_seq.substr(0, k_len), k_len, true);
    auto iter1 = code2kmer.find(kmer_code1);
    while (start_pos1 + k_len < seq_size)
    {
        if ((iter1 = code2kmer.find(kmer_code1)) != code2kmer.cend() ||
            (!stranded && (iter1 = code2kmer.find(GetRC(kmer_code1, k_len))) != code2kmer.cend()))
        {
            break;
        }
        kmer_code1 = NextCode(kmer_code1, k_len, contig_seq[start_pos1 + k_len]);
        start_pos1++;
    }
    if (start_pos1 + k_len >= seq_size) // not able to find a left k-mer to start
    {
        max_kmer1 = max_kmer2 = "NONE";
        return kMaxDistance;
    }
    no_norm ? iter1->second.GetCountVect(left_count_vect, idx_file, nb_count) : iter1->second.GetCountVect(left_count_vect, idx_file, nb_count, nf_vect);

    // Traverse the whole sequence //
    float max_dist = kMinDistance, dist_x;
    size_t start_pos2 = start_pos1;
    uint64_t kmer_code2 = kmer_code1;
    auto iter2 = iter1;
    while (start_pos1 + k_len < seq_size)
    {
        do
        {
            kmer_code2 = NextCode(kmer_code2, k_len, contig_seq[start_pos2 + k_len]);
            start_pos2++;
            if ((iter2 = code2kmer.find(kmer_code2)) != code2kmer.cend() ||
                (!stranded && (iter2 = code2kmer.find(GetRC(kmer_code2, k_len))) != code2kmer.cend()))
            {
                break;
            }
        } while (start_pos2 + k_len <= seq_size && start_pos2 - start_pos1 <= max_shift);
        if (start_pos2 + k_len > seq_size || start_pos2 - start_pos1 > max_shift) // not able to find a right k-mer within allowed zone
        {
            Int2Seq(max_kmer1, kmer_code1, k_len);
            max_kmer2 = "NONE";
            return kMaxDistance;
        }
        no_norm ? iter2->second.GetCountVect(right_count_vect, idx_file, nb_count) : iter2->second.GetCountVect(right_count_vect, idx_file, nb_count, nf_vect);
        if ((eval_method == "pearson" && (dist_x = CalcPearsonDist(left_count_vect, right_count_vect)) > max_dist) ||
            (eval_method == "spearman" && (dist_x = CalcSpearmanDist(left_count_vect, right_count_vect)) > max_dist) ||
            (eval_method == "mac" && (dist_x = CalcMACDist(left_count_vect, right_count_vect)) > max_dist)) // short-circuit operator
        {
            max_dist = dist_x;
            max_kmer_code1 = kmer_code1;
            max_kmer_code2 = kmer_code2;
        }
        start_pos1 = start_pos2;
        kmer_code1 = kmer_code2;
        iter1 = iter2;
        left_count_vect.swap(right_count_vect); // same as left_count_vect = right_count_vect, but in constant complexity
    }

    Int2Seq(max_kmer1, max_kmer_code1, k_len);
    Int2Seq(max_kmer2, max_kmer_code2, k_len);
    return max_dist;
}

const std::vector<float> &EvaluateSeqCount(std::vector<float> &count_vect, const std::string &contig_seq, const std::vector<double> &nf_vect,
                                           const code2kmer_t &code2kmer, std::ifstream &idx_file, const bool no_norm,
                                           const size_t k_len, const bool stranded, const std::string &eval_method)
{
    const size_t nb_count = count_vect.size(), seq_len = contig_seq.size();
    auto iter = code2kmer.cend();
    if (eval_method == "mean")
    {
        static std::vector<float> count_vect_x;
        size_t nb_mem_kmer = 0, kmer_code;
        kmer_code = Seq2Int(contig_seq.substr(0, k_len), k_len, true);
        for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
        {
            if (((iter = code2kmer.find(kmer_code)) != code2kmer.cend()) ||
                (!stranded && (iter = code2kmer.find(GetRC(kmer_code, k_len))) != code2kmer.cend()))
            {
                if (no_norm)
                {
                    iter->second.GetCountVect(count_vect_x, idx_file, nb_count);
                }
                else
                {
                    iter->second.GetCountVect(count_vect_x, idx_file, nb_count, nf_vect);
                }
                for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
                {
                    count_vect[i_smp] += count_vect_x[i_smp];
                }
                nb_mem_kmer++;
            }
            if (start_pos + k_len < seq_len)
            {
                kmer_code = NextCode(kmer_code, k_len, contig_seq[start_pos + k_len]);
            }
        }
        for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
        {
            count_vect[i_smp] /= nb_mem_kmer;
        }
    }
    else if (eval_method == "median")
    {
        static std::vector<std::vector<float>> kmer_count_mat;
        static std::vector<float> count_vect_x;
        kmer_count_mat.clear();
        size_t kmer_code = Seq2Int(contig_seq.substr(0, k_len), k_len, true);
        for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
        {
            if (((iter = code2kmer.find(kmer_code)) != code2kmer.cend()) ||
                (!stranded && (iter = code2kmer.find(GetRC(kmer_code, k_len))) != code2kmer.cend()))
            {
                kmer_count_mat.emplace_back(std::vector<float>(nb_count, 0));
                if (no_norm)
                {
                    iter->second.GetCountVect(kmer_count_mat.back(), idx_file, nb_count);
                }
                else
                {
                    iter->second.GetCountVect(kmer_count_mat.back(), idx_file, nb_count, nf_vect);
                }
            }
            if (start_pos + k_len < seq_len)
            {
                kmer_code = NextCode(kmer_code, k_len, contig_seq[start_pos + k_len]);
            }
        }
        const size_t nb_mem_kmer = kmer_count_mat.size();
        count_vect_x.resize(nb_mem_kmer);
        for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
        {
            for (size_t i_kmer(0); i_kmer < nb_mem_kmer; ++i_kmer)
            {
                count_vect_x[i_kmer] = kmer_count_mat[i_kmer][i_smp];
            }
            count_vect[i_smp] = CalcVectMedian(count_vect_x);
        }
    }
    return count_vect;
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string contig_fasta_path, eval_method, colname_list_path, idx_path, kmer_count_path, out_path;
    size_t k_len, max_shift = std::numeric_limits<size_t>::max();
    bool stranded(true), no_norm(false), contig_name(false);

    ParseOptions(argc, argv, k_len, contig_fasta_path, eval_method, idx_path, stranded, colname_list_path, no_norm, contig_name, out_path, kmer_count_path);
    PrintRunInfo(k_len, contig_fasta_path, eval_method, idx_path, stranded, colname_list_path, no_norm, contig_name, out_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    TabHeader tab_header(colname_list_path);
    code2kmer_t code2kmer;
    std::vector<double> nf_vect;
    ScanCountTable(tab_header, code2kmer, nf_vect, k_len, stranded, kmer_count_path, idx_path);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    fastaVect_t fasta_vect;
    EstablishSeqListFromMultilineFasta(fasta_vect, contig_fasta_path);
    std::cerr << "Number of sequence for evaluation: " << fasta_vect.size() << std::endl;

    if (eval_method == "pearson" || eval_method == "spearman" || eval_method == "mac")
    {
        std::cout << "contig\teval_dist\tkmer1\tkmer2" << std::endl;
    }
    else if (eval_method == "mean" || eval_method == "median")
    {
        std::cout << "contig";
        for (size_t i(0); i < tab_header.GetNbCol(); ++i)
        {
            if (tab_header.IsColCount(i))
            {
                std::cout << "\t" << tab_header.GetColNameAt(i);
            }
        }
        std::cout << std::endl;
    }
    std::string contig_seq, max_kmer1, max_kmer2;
    std::ifstream idx_file(idx_path);
    if (!idx_file.is_open())
    {
        throw std::domain_error("could not open index file " + idx_path);
    }
    for (const auto &seq_elem_ptr : fasta_vect)
    {
        contig_seq = seq_elem_ptr->second;
        if (contig_name)
        {
            std::cout << seq_elem_ptr->first;
        }
        else
        {
            std::cout << contig_seq;
        }
        if (eval_method == "pearson" || eval_method == "spearman" || eval_method == "mac")
        {
            std::cout << "\t" << EvaluateSeqDist(max_kmer1, max_kmer2, contig_seq, nf_vect, no_norm, k_len, stranded, eval_method, code2kmer, idx_file, tab_header.GetNbCount(), max_shift)
                      << "\t" << max_kmer1 << "\t" << max_kmer2 << std::endl;
        }
        else if (eval_method == "mean" || eval_method == "median")
        {
            static std::vector<float> count_vect(tab_header.GetNbCount(), 0);
            for (const float c : EvaluateSeqCount(count_vect, contig_seq, nf_vect, code2kmer, idx_file, no_norm, k_len, stranded, eval_method))
            {
                std::cout << "\t" << c;
            }
            std::cout << std::endl;
        }
    }
    std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;

    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}
