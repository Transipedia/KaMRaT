#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <limits>
#include <ctime>

#include "runinfo_files/query_runinfo.hpp"

const float kMinDistance = 0, kMaxDistance = 1;

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp
uint64_t GetRC(const uint64_t code, size_t k_length);                                 // in utils/seq_coding.cpp
uint64_t NextCode(uint64_t code, const size_t k_length, const char new_nuc);          // in utils/seq_coding.cpp

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, const std::string &idx_meta_path);   // in utils/index_loading.cpp
void LoadCodePosMap(std::map<uint64_t, size_t> &code_pos_map, const std::string &idx_pos_path); // in utils/index_loading.cpp

const std::vector<float> &GetMeanCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                           const std::vector<size_t> &mem_pos_vect); // in utils/index_loading.cpp
const std::vector<float> &GetMedianCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                             const std::vector<size_t> &mem_pos_vect); // in utils/index_loading.cpp

void PrintHeader(const std::vector<std::string> &colname_vect)
{
    std::cout << colname_vect[0];
    for (size_t i(1); i < colname_vect.size(); ++i)
    {
        std::cout << "\t" << colname_vect[i];
    }
    std::cout << std::endl;
}

void PrintQueryRes(const std::string &seq, const std::string &query_mthd,
                   std::ifstream &idx_mat, const std::map<uint64_t, size_t> &code_pos_map,
                   const size_t nb_smp, const size_t k_len, const bool stranded, const bool with_absent)
{
    static std::vector<size_t> mem_pos_vect;
    static std::vector<float> count_vect;
    const size_t seq_len = seq.size();
    mem_pos_vect.reserve(seq_len - k_len + 1);
    std::map<uint64_t, size_t>::const_iterator iter;
    uint64_t kmer_code(Seq2Int(seq.substr(0, k_len), k_len, true));
    for (size_t start_pos(0); start_pos + k_len - 1 < seq_len; ++start_pos)
    {
        if (((iter = code_pos_map.find(kmer_code)) != code_pos_map.cend()) ||
            (!stranded && (iter = code_pos_map.find(GetRC(kmer_code, k_len))) != code_pos_map.cend()))
        {
            mem_pos_vect.emplace_back(iter->second);
        }
        if (start_pos + k_len < seq_len)
        {
            kmer_code = NextCode(kmer_code, k_len, seq[start_pos + k_len]);
        }
    }
    if (mem_pos_vect.empty() && with_absent)
    {
        std::cout << seq;
        for (size_t i(0); i < nb_smp; ++i)
        {
            std::cout << "\t0";
        }
        std::cout << std::endl;
    }
    else if (!mem_pos_vect.empty())
    {
        if (query_mthd == "mean")
        {
            GetMeanCountVect(count_vect, idx_mat, nb_smp, mem_pos_vect);
        }
        else
        {
            GetMedianCountVect(count_vect, idx_mat, nb_smp, mem_pos_vect);
        }
        std::cout << seq;
        for (const float x : count_vect)
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
    }
    mem_pos_vect.clear();
    count_vect.clear();
}

int QueryMain(int argc, char **argv)
{
    QueryWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, seq_file_path, query_mtd, out_path;
    size_t nb_smp, k_len;
    bool stranded(true), with_absent(false);
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, seq_file_path, query_mtd, with_absent, out_path);
    PrintRunInfo(idx_dir, seq_file_path, query_mtd, with_absent, out_path);

    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, idx_dir + "/idx-meta.bin");
    if (k_len == 0)
    {
        throw std::invalid_argument("KaMRaT-query relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }
    std::map<uint64_t, size_t> code_pos_map;
    LoadCodePosMap(code_pos_map, idx_dir + "/idx-pos.bin");
    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

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
    PrintHeader(colname_vect);

    std::ifstream seq_list_file(seq_file_path);
    if (!seq_list_file.is_open())
    {
        throw std::invalid_argument("cannot open file: " + seq_file_path);
    }
    std::string seq, line;
    bool is_first_line(true);
    size_t nb_seq(0);
    while (std::getline(seq_list_file, line))
    {
        if (line[0] == '>') // if a header line (i.e. a new sequence)
        {
            if (!is_first_line)
            {
                PrintQueryRes(seq, query_mtd, idx_mat, code_pos_map, nb_smp, k_len, stranded, with_absent);
                seq.clear(); // prepare for the next sequence
                nb_seq++;
            }
            is_first_line = false;
        }
        else
        {
            seq += line;
        }
    }
    PrintQueryRes(seq, query_mtd, idx_mat, code_pos_map, nb_smp, k_len, stranded, with_absent); // query the last sequence
    seq_list_file.close();
    std::cerr << "Number of sequence for evaluation: " << nb_seq << std::endl;

    idx_mat.close();
    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }

    std::cerr << "Query and print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}
