#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <limits>
#include <ctime>

#include "query/query_runinfo.hpp"

const float kMinDistance = 0, kMaxDistance = 1;

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp
uint64_t GetRC(const uint64_t code, size_t k_length);                                 // in utils/seq_coding.cpp
uint64_t NextCode(uint64_t code, const size_t k_length, const char new_nuc);          // in utils/seq_coding.cpp

const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y);  // in utils/vect_opera.cpp
const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y); // in utils/vect_opera.cpp
const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y);      // in utils/vect_opera.cpp

double QueryDist(std::string &kmer1, std::string &kmer2, std::ifstream &idx_mat,
                 const size_t nb_smp, const size_t k_len, const bool stranded, const size_t max_shift,
                 const std::unordered_map<uint64_t, size_t> &code_pos_map, const std::string &seq, const std::string &query_mtd)
{
    if (seq.size() == k_len && code_pos_map.find(Seq2Int(seq, k_len, stranded)) != code_pos_map.cend()) // if seq is a k-mer in k-mer count table
    {
        kmer1 = kmer2 = seq;
        return kMinDistance;
    }
    else if (seq.size() <= k_len) // if seq is a k-mer but not in k-mer count table or if it's shorter than a k-mer
    {
        kmer1 = kmer2 = "NONE";
        return kMaxDistance;
    }

    static std::vector<float> left_count_vect, right_count_vect;
    size_t seq_size = seq.size(), kmer_code1, kmer_code2;

    // Start condition: left k-mer as the first findable one in the k-mer count table //
    size_t start_pos1 = 0;
    uint64_t kmer_code1 = Seq2Int(seq.substr(0, k_len), k_len, true);
    auto iter1 = code_pos_map.find(kmer_code1);
    while (start_pos1 + k_len < seq_size)
    {
        if ((iter1 = code_pos_map.find(kmer_code1)) != code_pos_map.cend() ||
            (!stranded && (iter1 = code_pos_map.find(GetRC(kmer_code1, k_len))) != code_pos_map.cend()))
        {
            break;
        }
        kmer_code1 = NextCode(kmer_code1, k_len, seq[start_pos1 + k_len]);
        start_pos1++;
    }
    if (start_pos1 + k_len >= seq_size) // not able to find a left k-mer to start
    {
        kmer1 = kmer2 = "NONE";
        return kMaxDistance;
    }
    GetCountVect(left_count_vect, idx_mat, iter1->second, nb_smp);
    idx_mat >> kmer1;

    // Traverse the whole sequence //
    float max_dist = kMinDistance, dist_x;
    size_t start_pos2 = start_pos1;
    uint64_t kmer_code2 = kmer_code1;
    auto iter2 = iter1;
    while (start_pos1 + k_len < seq_size)
    {
        do
        {
            kmer_code2 = NextCode(kmer_code2, k_len, seq[start_pos2 + k_len]);
            start_pos2++;
            if ((iter2 = code_pos_map.find(kmer_code2)) != code_pos_map.cend() ||
                (!stranded && (iter2 = code_pos_map.find(GetRC(kmer_code2, k_len))) != code_pos_map.cend()))
            {
                break;
            }
        } while (start_pos2 + k_len <= seq_size && start_pos2 - start_pos1 <= max_shift);
        if (start_pos2 + k_len > seq_size || start_pos2 - start_pos1 > max_shift) // not able to find a right k-mer within allowed zone
        {
            kmer2 = "NONE";
            return kMaxDistance;
        }
        GetCountVect(right_count_vect, idx_mat, iter2->second, nb_smp);
        if ((query_mtd == "pearson" && (dist_x = CalcPearsonDist(left_count_vect, right_count_vect)) > max_dist) ||
            (query_mtd == "spearman" && (dist_x = CalcSpearmanDist(left_count_vect, right_count_vect)) > max_dist) ||
            (query_mtd == "mac" && (dist_x = CalcMACDist(left_count_vect, right_count_vect)) > max_dist)) // short-circuit operator
        {
            max_dist = dist_x;
            kmer_code1 = kmer_code1;
            kmer_code2 = kmer_code2;
        }
        start_pos1 = start_pos2;
        kmer_code1 = kmer_code2;
        iter1 = iter2;
        left_count_vect.swap(right_count_vect); // same as left_count_vect = right_count_vect, but in constant complexity
    }

    Int2Seq(kmer1, kmer_code1, k_len);
    Int2Seq(kmer2, kmer_code2, k_len);
    return max_dist;
}

void ScanQuery(const std::string &seq_file_path, std::ifstream &idx_mat, const std::unordered_map<uint64_t, size_t> &code_pos_map,
               const std::string &query_mtd, const bool out_name)
{
    std::ifstream contig_list_file(seq_file_path);
    if (!contig_list_file.is_open())
    {
        throw std::invalid_argument("cannot open file: " + seq_file_path);
    }
    std::string seq_name, seq, line, kmer1, kmer2;
    bool is_first_line(true);
    while (std::getline(contig_list_file, line))
    {
        if (line[0] == '>') // if a header line (i.e. a new sequence)
        {
            if (!is_first_line)
            {
                if (query_mtd == "pearson" || query_mtd == "spearman" || query_mtd == "mac")
                {
                    QueryDist(kmer1, kmer2, idx_mat, code_pos_map, seq, query_mtd);
                }
                seq.clear(); // prepare for the next sequence
            }
            PrintDist(out_name ? seq_name : seq, kmer1, kmer2);
            seq_name = line.substr(1); // record the next sequence name
            is_first_line = false;
        }
        else
        {
            seq += line;
        }
    }
    seq_vect.emplace_back(std::make_pair(seq_name, seq)); // cache the last sequence
    contig_list_file.close();
}

// const float EvaluateSeqDist(std::string &kmer1, std::string &kmer2,
//                             const std::string &seq, const std::vector<double> &nf_vect, const bool no_norm,
//                             const size_t k_len, const bool stranded, const std::string &query_mtd,
//                             const code2kmer_t &code_pos_map, std::ifstream &idx_file, const size_t nb_count, const size_t max_shift)
// {
//
// }

// const std::vector<float> &EvaluateSeqCount(std::vector<float> &count_vect, const std::string &contig_seq, const std::vector<double> &nf_vect,
//                                            const code2kmer_t &code2kmer, std::ifstream &idx_file, const bool no_norm,
//                                            const size_t k_len, const bool stranded, const std::string &query_mtd)
// {
//     const size_t nb_count = count_vect.size(), seq_len = contig_seq.size();
//     auto iter = code2kmer.cend();
//     if (query_mtd == "mean")
//     {
//         static std::vector<float> count_vect_x;
//         size_t nb_mem_kmer = 0, kmer_code;
//         kmer_code = Seq2Int(contig_seq.substr(0, k_len), k_len, true);
//         for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
//         {
//             if (((iter = code2kmer.find(kmer_code)) != code2kmer.cend()) ||
//                 (!stranded && (iter = code2kmer.find(GetRC(kmer_code, k_len))) != code2kmer.cend()))
//             {
//                 if (no_norm)
//                 {
//                     iter->second.GetCountVect(count_vect_x, idx_file, nb_count);
//                 }
//                 else
//                 {
//                     iter->second.GetCountVect(count_vect_x, idx_file, nb_count, nf_vect);
//                 }
//                 for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
//                 {
//                     count_vect[i_smp] += count_vect_x[i_smp];
//                 }
//                 nb_mem_kmer++;
//             }
//             if (start_pos + k_len < seq_len)
//             {
//                 kmer_code = NextCode(kmer_code, k_len, contig_seq[start_pos + k_len]);
//             }
//         }
//         for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
//         {
//             count_vect[i_smp] /= nb_mem_kmer;
//         }
//     }
//     else if (query_mtd == "median")
//     {
//         static std::vector<std::vector<float>> kmer_count_mat;
//         static std::vector<float> count_vect_x;
//         kmer_count_mat.clear();
//         size_t kmer_code = Seq2Int(contig_seq.substr(0, k_len), k_len, true);
//         for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
//         {
//             if (((iter = code2kmer.find(kmer_code)) != code2kmer.cend()) ||
//                 (!stranded && (iter = code2kmer.find(GetRC(kmer_code, k_len))) != code2kmer.cend()))
//             {
//                 kmer_count_mat.emplace_back(std::vector<float>(nb_count, 0));
//                 if (no_norm)
//                 {
//                     iter->second.GetCountVect(kmer_count_mat.back(), idx_file, nb_count);
//                 }
//                 else
//                 {
//                     iter->second.GetCountVect(kmer_count_mat.back(), idx_file, nb_count, nf_vect);
//                 }
//             }
//             if (start_pos + k_len < seq_len)
//             {
//                 kmer_code = NextCode(kmer_code, k_len, contig_seq[start_pos + k_len]);
//             }
//         }
//         const size_t nb_mem_kmer = kmer_count_mat.size();
//         if (nb_mem_kmer == 0)
//         {
//             count_vect.assign(nb_count, 0);
//         }
//         else
//         {
//             count_vect_x.resize(nb_mem_kmer);
//             for (size_t i_smp(0); i_smp < nb_count; ++i_smp)
//             {
//                 for (size_t i_kmer(0); i_kmer < nb_mem_kmer; ++i_kmer)
//                 {
//                     count_vect_x[i_kmer] = kmer_count_mat[i_kmer][i_smp];
//                 }
//                 count_vect[i_smp] = CalcVectMedian(count_vect_x);
//             }
//             count_vect_x.clear();
//         }
//     }
//     return count_vect;
// }

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded, std::vector<std::string> &colname_vect, const std::string &idx_meta_path); // in utils/index_loading.cpp
void LoadCodePosMap(std::unordered_map<uint64_t, size_t> &code_pos_map, const std::string &idx_pos_path);                                        // in utils/index_loading.cpp

int main(int argc, char **argv)
{
    QueryWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, seq_file_path, query_mtd, out_path;
    size_t nb_smp, k_len, max_shift = std::numeric_limits<size_t>::max();
    bool stranded(true), out_name(false);
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, seq_file_path, query_mtd, max_shift, out_name, out_path);
    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, idx_dir + "/idx-meta.bin");
    if (k_len == 0)
    {
        throw std::invalid_argument("KaMRaT-query relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }
    PrintRunInfo(idx_dir, seq_file_path, query_mtd, max_shift, out_name, out_path);
    std::unordered_map<uint64_t, size_t> code_pos_map;
    LoadCodePosMap(code_pos_map, idx_dir + "/idx-pos.bin");
    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    ScanQuery(seq_file_path);
    std::cerr << "Number of sequence for evaluation: " << fasta_vect.size() << std::endl;

    // std::ofstream out_file;
    // if (!out_path.empty())
    // {
    //     out_file.open(out_path);
    //     if (!out_file.is_open())
    //     {
    //         throw std::domain_error("cannot open file: " + out_path);
    //     }
    // }
    // auto backup_buf = std::cout.rdbuf();
    // if (!out_path.empty()) // output to file if a path is given, to screen if not
    // {
    //     std::cout.rdbuf(out_file.rdbuf());
    // }
    // if (query_mtd == "pearson" || query_mtd == "spearman" || query_mtd == "mac")
    // {
    //     std::cout << "contig\teval_dist\tkmer1\tkmer2" << std::endl;
    // }
    // else if (query_mtd == "mean" || query_mtd == "median")
    // {
    //     std::cout << "contig";
    //     for (size_t i(0); i < tab_header.GetNbCol(); ++i)
    //     {
    //         if (tab_header.IsColCount(i))
    //         {
    //             std::cout << "\t" << tab_header.GetColNameAt(i);
    //         }
    //     }
    //     std::cout << std::endl;
    // }
    // std::string contig_seq, max_kmer1, max_kmer2;
    // std::ifstream idx_file(idx_path);
    // if (!idx_file.is_open())
    // {
    //     throw std::domain_error("could not open index file " + idx_path);
    // }
    // for (const std::pair<std::string, std::string> &seq_elem : fasta_vect)
    // {
    //     contig_seq = seq_elem.second;
    //     if (out_name)
    //     {
    //         std::cout << seq_elem.first;
    //     }
    //     else
    //     {
    //         std::cout << contig_seq;
    //     }
    //     if (query_mtd == "pearson" || query_mtd == "spearman" || query_mtd == "mac")
    //     {
    //         std::cout << "\t" << EvaluateSeqDist(max_kmer1, max_kmer2, contig_seq, nf_vect, no_norm, k_len, stranded, query_mtd, code2kmer, idx_file, tab_header.GetNbCount(), max_shift)
    //                   << "\t" << max_kmer1 << "\t" << max_kmer2 << std::endl;
    //     }
    //     else if (query_mtd == "mean" || query_mtd == "median")
    //     {
    //         static std::vector<float> count_vect(tab_header.GetNbCount(), 0);
    //         for (const float c : EvaluateSeqCount(count_vect, contig_seq, nf_vect, code2kmer, idx_file, no_norm, k_len, stranded, query_mtd))
    //         {
    //             std::cout << "\t" << c;
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    // std::cout.rdbuf(backup_buf);
    // if (out_file.is_open())
    // {
    //     out_file.close();
    // }
    // std::cerr << "Contig evaluation finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;

    // std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    // return EXIT_SUCCESS;

    std::ofstream out_file;
    if (!out_path.empty())
    {
        out_file.open(out_path);
        if (!out_file.is_open())
        {
            throw std::domain_error("cannot open file: " + out_path);
        }
    }

    auto backup_buf = std::cout.rdbuf();
    if (!out_path.empty()) // output to file if a path is given, to screen if not
    {
        std::cout.rdbuf(out_file.rdbuf());
    }
    std::cout << colname_vect[0];
    for (size_t i_col(1); i_col <= nb_smp; ++i_col)
    {
        std::cout << "\t" << colname_vect[i_col];
    }
    std::cout << std::endl;

    idx_mat.close();
    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}
