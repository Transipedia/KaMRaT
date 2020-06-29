#include <iostream>
#include <string>
#include <unordered_set>
#include <fstream>
#include <sstream>

#include "utils/seq_coding.hpp"
#include "data_struct/seq_elem.hpp"
#include "run_info_parser/mask.hpp"

void MakeMask(std::unordered_set<uint64_t> &kmer_mask,
              const std::string &mask_file_path,
              const size_t k_length,
              const bool stranded)
{
    std::ifstream contig_list_file(mask_file_path);
    if (!contig_list_file.is_open())
    {
        throw std::domain_error("contig fasta file " + mask_file_path + "was not found");
    }
    std::string seq, line;
    for (std::getline(contig_list_file, line); std::getline(contig_list_file, line); true) // ignore the first header line
    {
        if (line[0] == '>' || contig_list_file.eof()) // reading the last line
        {
            for (size_t start_pos(0); start_pos + k_length <= seq.size(); ++start_pos)
            {
                std::string kmer = seq.substr(start_pos, k_length);
                kmer_mask.insert(Seq2Int(kmer, k_length, stranded));
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

void MaskOutput(std::string &&count_tab_path,
                const std::unordered_set<uint64_t> &kmer_mask,
                const size_t k_length,
                const bool stranded,
                const bool reverse_mask)
{
    std::ifstream kmer_count_file(count_tab_path);
    std::string line;
    if (!kmer_count_file.is_open())
    {
        std::cerr << "ERROR: Filter list file " << count_tab_path << " does not exist..." << std::endl;
        exit(EXIT_FAILURE);
    }
    // dealing with header line //
    std::getline(kmer_count_file, line);
    std::cout << line << std::endl;
    // dealing with other lines //
    for (std::getline(kmer_count_file, line); !kmer_count_file.eof(); std::getline(kmer_count_file, line))
    {
        std::istringstream conv(line);
        std::string kmer;
        conv >> kmer;
        bool is_in_mask = (kmer_mask.find(Seq2Int(kmer, k_length, stranded)) != kmer_mask.cend());
        if (is_in_mask && reverse_mask) // reverse_mask = to_select
        {
            std::cout << line << std::endl;
        }
        else if (!is_in_mask && !reverse_mask)
        {
            std::cout << line << std::endl;
        }
    }
    kmer_count_file.close();
}

int main(int argc, char **argv)
{
    std::string mask_file_path, count_tab_path;
    bool stranded(true), reverse_mask(false);
    size_t k_length;

    ParseOptions(argc, argv, mask_file_path, stranded, k_length, reverse_mask, count_tab_path);
    
    std::unordered_set<uint64_t> kmer_mask;
    MakeMask(kmer_mask, mask_file_path, k_length, stranded);
    MaskOutput(std::move(count_tab_path), kmer_mask, k_length, stranded, reverse_mask);
    return EXIT_SUCCESS;
}