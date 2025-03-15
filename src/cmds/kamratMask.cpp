#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <ctime>

#include "mask_runinfo.hpp"
#include "index_loading.hpp"
#include "seq_coding.hpp"

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

void ScanPrint(std::ifstream &idx_pos, std::ifstream &idx_mat,
               const std::unordered_set<uint64_t> &kmer2sel,
               const std::unordered_set<uint64_t> &kmer2sup,
               const bool with_counts, const size_t nb_smp,
               const std::string &value_mode)
{
    std::string kmer_seq;
    std::vector<float> count_vect;
    size_t code, pos;
    while (idx_pos.read(reinterpret_cast<char *>(&code), sizeof(size_t)) && idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        const bool is2sel = kmer2sel.empty() || (kmer2sel.find(code) != kmer2sel.cend());
        const bool is2sup = (!kmer2sup.empty()) && (kmer2sup.find(code) != kmer2sup.cend());
        if (is2sel && !is2sup)
        {
            GetCountVect(count_vect, idx_mat, pos, nb_smp);
            idx_mat >> kmer_seq;
            std::cout << kmer_seq;
            if (with_counts)
            {
                for (const float x : count_vect)
                {
                    if (value_mode == "float")
                    {
                        std::cout << "\t" << x;
                    }
                    else // value_mode == "int"
                    {
                        std::cout << "\t" << static_cast<size_t>(x + 0.5);
                    }
                }
            }
            else // output the binary file for other modules' `-with` argument
            {
                std::cout << "\t0\t1\t";
                std::cout.write(reinterpret_cast<char *>(&pos), sizeof(size_t));
            }
            std::cout << std::endl;
        }
    }
}

int MaskMain(int argc, char **argv)
{
    MaskWelcome();

    std::clock_t begin_time = clock();
    std::string idx_dir, seq2sel_path, seq2sup_path, out_fmt("tab"), out_path, value_mode("int");
    size_t nb_smp, k_len(0);
    bool stranded;
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, seq2sel_path, seq2sup_path, out_fmt, out_path, value_mode);
    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, idx_dir + "/idx-meta.bin");
    PrintRunInfo(idx_dir, k_len, stranded, seq2sel_path, seq2sup_path, out_fmt, out_path, value_mode);
    if (k_len == 0)
    {
        throw std::invalid_argument("KaMRaT-mask relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }

    std::unordered_set<uint64_t> kmer2sel, kmer2sup;
    if (!seq2sel_path.empty()) // index k-mers to be selected
    {
        MakeMask(kmer2sel, seq2sel_path, k_len, stranded);
    }
    if (!seq2sup_path.empty()) // index k-mers to be suppressed
    {
        MakeMask(kmer2sup, seq2sup_path, k_len, stranded);
    }

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
    if (out_fmt == "tab") // either `tab` or `bin`
    {
        std::cout << colname_vect[0];
        for (size_t i_col(1); i_col <= nb_smp; ++i_col)
        {
            std::cout << "\t" << colname_vect[i_col];
        }
        std::cout << std::endl;
    }
    ScanPrint(idx_pos, idx_mat, kmer2sel, kmer2sup, out_fmt == "tab", nb_smp, value_mode);
    idx_pos.close(), idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}