#include <iostream>
#include <string>
#include <unordered_set>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
    // std::getline(contig_list_file, line); // ignore the first header line
    while (true)
    {
        std::getline(contig_list_file, line);
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

    size_t pos = count_tab_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && count_tab_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    // dealing with header line //
    std::getline(kmer_count_instream, line);
    std::cout << line << std::endl;
    // dealing with other lines //
    for (std::getline(kmer_count_instream, line); !kmer_count_instream.eof(); std::getline(kmer_count_instream, line))
    {
        std::istringstream conv(line);
        std::string kmer;
        conv >> kmer;
        if (kmer.size() != k_length)
        {
            throw std::domain_error("k-mer length in table (" + std::to_string(kmer.size()) + ")" +
                                    " not equal to the given kmer length (" + std::to_string(k_length) + ")");
        }
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

int MaskMain(int argc, char **argv)
{
    std::string mask_file_path, count_tab_path;
    bool stranded(true), reverse_mask(false);
    size_t k_length(0);

    ParseOptions(argc, argv, k_length, mask_file_path, stranded, reverse_mask, count_tab_path);
    PrintRunInfo(k_length, mask_file_path, stranded, reverse_mask);

    std::unordered_set<uint64_t> kmer_mask;
    MakeMask(kmer_mask, mask_file_path, k_length, stranded);
    MaskOutput(std::move(count_tab_path), kmer_mask, k_length, stranded, reverse_mask);
    return EXIT_SUCCESS;
}