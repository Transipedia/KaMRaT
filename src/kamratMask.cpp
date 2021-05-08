#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <ctime>

#include "mask/mask_runinfo.hpp"

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path); // in utils/index_loading.cpp
const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat,
                                       const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp

void MakeMask(std::unordered_set<uint64_t> &kmer_mask, const std::string &mask_file_path, const size_t k_len, const bool stranded)
{
    std::ifstream contig_list_file(mask_file_path);
    if (!contig_list_file.is_open())
    {
        throw std::domain_error("contig fasta file " + mask_file_path + " was not found");
    }
    std::string seq, line;
    // std::getline(contig_list_file, line); // ignore the first header line
    while (true)
    {
        std::getline(contig_list_file, line);
        if (line[0] == '>' || contig_list_file.eof()) // reading the last line
        {
            for (size_t start_pos(0); start_pos + k_len <= seq.size(); ++start_pos)
            {
                std::string kmer = seq.substr(start_pos, k_len);
                kmer_mask.insert(Seq2Int(kmer, k_len, stranded));
            }
            seq.clear();
            if (contig_list_file.eof())
            {
                break;
            }
        }
        else
        {
            seq += line;
        }
    }
    contig_list_file.close();
}

void ScanPrint(std::ifstream &idx_pos, std::ifstream &idx_mat, const std::unordered_set<uint64_t> &kmer_mask,
               const bool reverse_mask, const bool with_counts, const size_t nb_smp)
{
    std::string kmer_seq;
    std::vector<float> count_vect;
    size_t code, pos;
    while (idx_pos.read(reinterpret_cast<char *>(&code), sizeof(size_t)) && idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        const bool is_in_mask = (kmer_mask.find(code) == kmer_mask.cend());
        if (is_in_mask == reverse_mask) // (is_in_mask && reverse_mask) || (!is_in_mask && !reverse_mask)
        {
            GetCountVect(count_vect, idx_mat, pos, nb_smp);
            idx_mat >> kmer_seq;
            std::cout << kmer_seq;
            if (with_counts)
            {
                for (const float x : count_vect)
                {
                    std::cout << "\t" << x;
                }
            }
            std::cout << std::endl;
        }
    }
}

int MaskMain(int argc, char **argv)
{
    MaskWelcome();

    std::clock_t begin_time = clock();
    std::string idx_dir, mask_file_path, out_path;
    size_t nb_smp, k_len;
    bool stranded, reverse_mask(false), with_counts(false);
    std::vector<std::string> colname_vect;
    std::vector<double> _smp_sum_vect; // _smp_sum_vect not needed
    ParseOptions(argc, argv, idx_dir, mask_file_path, reverse_mask, out_path, with_counts);
    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, _smp_sum_vect, idx_dir + "/idx-meta.bin");
    PrintRunInfo(idx_dir, k_len, stranded, mask_file_path, reverse_mask, out_path, with_counts);
    if (k_len == 0)
    {
        throw std::invalid_argument("KaMRaT-mask relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }

    std::unordered_set<uint64_t> kmer_mask;
    MakeMask(kmer_mask, mask_file_path, k_len, stranded);

    std::ifstream idx_pos(idx_dir + "/idx-pos.bin"), idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_pos.is_open() || !idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-pos or index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
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
    if (with_counts)
    {
        for (const auto &s : colname_vect)
        {
            std::cout << "\t" << s;
        }
        std::cout << std::endl;
    }
    ScanPrint(idx_pos, idx_mat, kmer_mask, reverse_mask, with_counts, nb_smp);
    idx_pos.close(), idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}