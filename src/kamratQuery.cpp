#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <numeric>
#include <ctime>

#include "utils/utils.hpp"
#include "utils/sample_info.hpp"
#include "utils/statistics.hpp"
#include "utils/seq_coding.hpp"
#include "query/parse_opt_print_info.hpp"
#include "query/contig_sample_info.hpp"

const size_t ScanCountTable(std::unordered_map<uint64_t, size_t> &kmer_disk_pos,
                            const std::string &kmer_count_path,
                            const SampleInfo &sample_info,
                            const std::string &count_index_path,
                            const bool stranded)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        std::cerr << "ERROR: k-mer count file " << kmer_count_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    // dealing with the header line for parsing column names //
    std::getline(kmer_count_file, line);
    std::vector<int> header_term_nat;
    size_t nb_sample = ParseHeader(header_term_nat, line, sample_info, "");
    // dealing with the count lines //
    std::ofstream count_index_file(count_index_path);
    for (size_t i_line(1); std::getline(kmer_count_file, line); ++i_line)
    {
        std::istringstream conv(line);
        std::string term, seq;
        std::vector<size_t> sample_counts;
        float rep_value;
        conv >> seq; // first column is sequence
        for (size_t i(0); conv >> term; ++i)
        {
            if (header_term_nat.at(i) >= 0)
            {
                sample_counts.push_back(static_cast<size_t>(std::stod(term) + 0.5)); // convert string to double and round it to integer
            }
        }
        if (!kmer_disk_pos.insert({Seq2Int(seq, seq.size(), stranded), WriteCountToIndex(count_index_file, sample_counts)}).second)
        {
            std::cerr << "ERROR: found duplicated k-mer sequence " << seq << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    count_index_file.close();
    kmer_count_file.close();

    return nb_sample;
}

void PrintQueryCount(std::ofstream &out_file,//////////////////////////////////// 问题：无法处理eliminate操作！！！！！！！！！！！！！！
                     std::ifstream &index_file,
                     const std::string &contig_seq,
                     const size_t nb_sample,
                     const size_t k_length,
                     const bool stranded,
                     const kmer_disk_pos_t kmer_pos_dict)
{
    std::vector<std::vector<size_t>> sample_kmer_count(nb_sample);
    size_t seq_len(contig_seq.size()), nb_kmers(1);
    std::string kmer_seq(contig_seq.substr(0, k_length));
    uint64_t kmer_code = Seq2Int(kmer_seq, k_length, stranded);
    auto iter = kmer_pos_dict.find(kmer_code);
    // get other k-mers one by one //
    std::vector<size_t> kmer_counts;
    if (iter != kmer_pos_dict.cend())
    {
        LoadCountFromIndex(index_file, kmer_counts, iter->second, nb_sample);
        for (size_t i(0); i < nb_sample; ++i)
        {
            sample_kmer_count.at(i).push_back(kmer_counts.at(i));
        }
    }
    for (size_t new_nuc_loc(k_length); new_nuc_loc < seq_len; ++new_nuc_loc)
    {
        kmer_code = next_kmer(k_length, contig_seq.at(new_nuc_loc), kmer_code, stranded);
        iter = kmer_pos_dict.find(kmer_code);
        if (iter != kmer_pos_dict.cend())
        {
            std::vector<size_t> count_x;
            LoadCountFromIndex(index_file, count_x, iter->second, nb_sample);
            for (size_t i(0); i < nb_sample; ++i)
            {
                sample_kmer_count.at(i).push_back(kmer_counts.at(i));
            }
        }
    }
}

void ScanContigList(const std::string &contig_list_path)
{
    std::ifstream contig_file(contig_list_path);
    if (!contig_file.is_open())
    {
        std::cerr << "ERROR: contig file " << contig_list_path << " was not found..." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(contig_file, line))
    {
        if (line[0] == '>')
        {
            continue;
        }
        std::string contig_seq = line.substr(0, line.find_first_of("\t \n"));
    }
    contig_file.close();
}

// void AddSampleCounts(ContigSampleInfo &contig_sample_info,
//                      std::ifstream &in_fs,
//                      const int k_length,
//                      const kmer_disk_pos_t &kmer_disk_pos,
//                      const std::vector<std::string> &col_name_vect,
//                      const std::vector<bool> &is_col_sample)
// {
//     counts.resize(nb_sample);
//     // get first k-mer //
//     size_t seq_len(seq.size()), nb_kmers(1);
//     std::string kmer_seq(seq.substr(0, k_length));
//     uint64_t kmer_code = Seq2Int(kmer_seq, k_length, stranded);
//     auto iter = kmer_pos_dict.find(kmer_code);
//     if (iter == kmer_pos_dict.cend())
//     {
//         std::cerr << "ERROR: k-mer not found in count index " << kmer_seq << std::endl;
//         exit(EXIT_FAILURE);
//     }
//     // get other k-mers one by one //
//     std::vector<size_t> sum;
//     LoadCountFromIndex(index_file, sum, iter->second, nb_sample);
//     for (size_t new_nuc_loc(k_length); new_nuc_loc < seq_len; ++new_nuc_loc)
//     {
//         kmer_code = next_kmer(k_length, seq.at(new_nuc_loc), kmer_code, stranded);
//         iter = kmer_pos_dict.find(kmer_code);
//         if (iter != kmer_pos_dict.cend())
//         {
//             std::vector<size_t> count_x;
//             LoadCountFromIndex(index_file, count_x, iter->second, nb_sample);
//             for (size_t i(0); i < nb_sample; ++i)
//             {
//                 sum.at(i) += count_x.at(i);
//             }
//             ++nb_kmers;
//         }
//     }
//     // turn sum to average //
//     for (size_t i = 0; i < nb_sample; ++i)
//     {
//         counts.at(i) = static_cast<size_t>(1.0 * sum.at(i) / nb_kmers + 0.5);
//     }

//     std::string contig_seq(contig_sample_info.GetContigSeq());
//     int max_nb_kmer(contig_seq.size() - k_length + 1);
//     for (int start_pos(0); start_pos < max_nb_kmer; ++start_pos)
//     {
//         std::string kmer_seq(contig_seq.substr(start_pos, k_length)), count_line;
//         auto iter = kmer_disk_pos.find(kmer_seq);
//         if (iter != kmer_disk_pos.cend())
//         {
//             GetStringLineFromDisk(count_line, in_fs, iter->second);
//             contig_sample_info.AddKmerCounts(count_line, col_name_vect, is_col_sample);
//         }
//     }
// }

int main(int argc, char *argv[])
{
    int k_length(31);
    std::string contig_list_path, colname_list_path, query_mode("extract"), estimate_method, kmer_count_path;
    bool stranded(true);
    ParseOptions(argc, argv, contig_list_path, k_length, stranded, colname_list_path, query_mode, estimate_method, kmer_count_path);
    PrintRunInfo(contig_list_path, k_length, stranded, colname_list_path, query_mode, estimate_method, kmer_count_path);

    SampleInfo sample_info;
    if (!colname_list_path.empty())
    {
        LoadSampleInfo(sample_info, colname_list_path);
    }

    kmer_disk_pos_t kmer_disk_pos;
    std::vector<std::string> col_name_vect;
    ScanCountTable(kmer_disk_pos, col_name_vect, kmer_count_path, stranded);

    std::vector<ContigSampleInfo> contig_sample_info_vect;
    MakeContigs(contig_sample_info_vect, contig_path);

    std::vector<bool> is_col_sample(col_name_vect.size(), true);
    is_col_sample[0] = false;
    if (!sample_name_set.empty())
    {
        for (int i(1); i < col_name_vect.size(); ++i)
        {
            is_col_sample[i] = (sample_name_set.find(col_name_vect[i]) != sample_name_set.cend());
        }
    }

    std::ifstream kmer_count_file(kmer_count_path);
    std::cout << col_name_vect[0];
    for (int i(1); i < col_name_vect.size(); ++i)
    {
        if (is_col_sample[i])
        {
            std::cout << "\t" << col_name_vect[i];
        }
    }
    std::cout << std::endl;
    for (auto &csi : contig_sample_info_vect)
    {
        AddSampleCounts(csi, kmer_count_file, k_length, kmer_disk_pos, col_name_vect, is_col_sample);
        csi.EstimatePrint(std::cout, &CalcVectMean, col_name_vect, is_col_sample);
    }
    kmer_count_file.close();

    return EXIT_SUCCESS;
}
