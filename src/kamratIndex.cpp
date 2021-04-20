#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "index/index_runinfo.hpp"

#define RESET "\033[0m"
#define BOLDYELLOW "\033[1m\033[33m"

/* ----------------------------------------------------------------- *\ 
 * idx-meta:                                                         *
 *   - sample number, k length (0 if general feature), strandedness  * 
 *   - header row indicating column names                            *
 *   - sample sum vector (binarized double vector)                   *
 * idx-ftpos:                                                        *
 *   - {feature string, position} ordered by feature                 *
 * idx-mat:                                                          *
 *   - feature counts (binarized float vector)                       *
\* ----------------------------------------------------------------- */

const uint64_t Seq2Int(const std::string &seq, const size_t k_len, const bool stranded); // seq_coding.cpp

const size_t CountColumn(std::ofstream &idx_meta, const std::string &line_str)
{
    size_t nb_smp(0);
    std::istringstream conv(line_str);
    std::string term;
    for (conv >> term; conv >> term; ++nb_smp) // count sample number, skipping the first column
    {
    }
    if (nb_smp == 0)
    {
        throw std::domain_error("input table parsing failed: sample number equals to 0");
    }
    return nb_smp;
}

const void IndexCount(std::ofstream &idx_mat, std::vector<double> &sum_vect, std::map<uint64_t, size_t> &code_pos_map,
                      const std::string &line_str, const size_t k_len, const bool stranded, const size_t nb_smp)
{
    static std::istringstream conv(line_str);
    static std::vector<float> count_vect;
    static std::string ft_name, term;
    static size_t ft_code(0), ft_pos;

    conv.str(line_str);
    ft_pos = static_cast<size_t>(idx_mat.tellp());
    for (conv >> ft_name; conv >> term; count_vect.push_back(std::stof(term))) // parse feature name and following count columns
    {
    }
    if (count_vect.size() != nb_smp) // check if all rows have same number of columns as the header row
    {
        throw std::length_error("sample numbers are not consistent: " + std::to_string(nb_smp) + " vs " + std::to_string(count_vect.size()));
    }
    for (size_t i_smp(0); i_smp < nb_smp; ++i_smp) // add count vectors together for eventual normalization
    {
        sum_vect[i_smp] += count_vect[i_smp];
    }
    if (k_len != 0) // if index in k-mer mode => ft_code calculated by Seq2Int
    {
        if (k_len != ft_name.size()) // check k-mer length
        {
            throw std::length_error("feature length checking failed: length of " + ft_name + " not equal to " + std::to_string(k_len));
        }
        ft_code = Seq2Int(ft_name, k_len, stranded);
    }
    else // if index in general mode => ft_code is serial number
    {
        ft_code++;
    }
    if (!code_pos_map.insert({ft_code, ft_pos}).second)
    {
        throw std::domain_error("unicity checking failed, an equivalent key already existed for feature: " + ft_name);
    }

    idx_mat.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(float)); // [idx_mat] feature count vector
    idx_mat << ft_name << std::endl;

    ft_name.clear();
    term.clear();
    count_vect.clear();
    conv.clear();
}

void IndexPos(std::ofstream &idx_pos, const std::map<uint64_t, size_t> &code_pos_map)
{
    uint64_t code;
    size_t pos;
    for (const auto &elem : code_pos_map) // [idx_pos] feature code and feature position, ordered by code
    {
        code = elem.first;
        pos = elem.second;
        idx_pos.write(reinterpret_cast<char *>(&code), sizeof(uint64_t));
        idx_pos.write(reinterpret_cast<char *>(&pos), sizeof(size_t));
    }
}

void ScanIndex(std::ofstream &idx_meta, std::ofstream &idx_pos, std::ofstream &idx_mat,
               std::istream &kmer_count_instream, const size_t k_len, const bool stranded)
{
    std::string line;
    std::getline(kmer_count_instream, line); // read header row in table
    size_t nb_smp = CountColumn(idx_meta, line);

    idx_meta << nb_smp << "\t" << k_len; // [idx_meta 1] sample number, k-mer length, and strandedness if applicable
    if (k_len != 0)
    {
        idx_meta << "\t" << (stranded ? 'T' : 'F');
    }
    idx_meta << std::endl;

    idx_meta << line << std::endl; // [idx_meta 2] the header row

    std::vector<double> sum_vect(nb_smp, 0);
    std::map<uint64_t, size_t> code_pos_map;
    while (std::getline(kmer_count_instream, line))
    {
        IndexCount(idx_mat, sum_vect, code_pos_map, line, k_len, stranded, nb_smp); // [idx_pos, idx_mat] (inside)
    }
    idx_meta << sum_vect[0]; // [idx_meta 3] sample sum vector
    for (size_t i(1); i < nb_smp; ++i)
    {
        idx_meta << "\t" << sum_vect[i];
    }
    idx_meta << std::endl;
    IndexPos(idx_pos, code_pos_map);
}

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path);                             // in utils/index_loading.cpp
void LoadPosVect(std::map<uint64_t, size_t> &code_pos_map, const std::string &idx_pos_path); // in utils/index_loading.cpp
const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp
void TestIndex(const std::string &idx_meta_path, const std::string &idx_pos_path, const std::string &idx_mat_path)
{
    std::ifstream idx_mat(idx_mat_path);
    std::string term;
    size_t nb_smp_all, k_len;
    bool stranded;
    std::vector<std::string> colnames_vect;
    std::vector<double> smp_sum_vect;
    LoadIndexMeta(nb_smp_all, k_len, stranded, colnames_vect, smp_sum_vect, idx_meta_path);
    std::cout << nb_smp_all << "\t" << k_len;
    if (k_len > 0)
    {
        std::cout << "\t" << (stranded ? "T" : "F");
    }
    std::cout << std::endl;
    for (double s : smp_sum_vect)
    {
        std::cout << s << "\t";
    }
    std::cout << std::endl;
    for (auto const &term : colnames_vect)
    {
        std::cout << term << "\t";
    }
    std::cout << std::endl;

    std::map<uint64_t, size_t> code_pos_map;
    LoadPosVect(code_pos_map, idx_pos_path);

    std::vector<float> count_vect;
    for (const auto &elem : code_pos_map)
    {
        GetCountVect(count_vect, idx_mat, elem.second, nb_smp_all);
        idx_mat >> term;
        std::cout << term;
        for (const auto x : count_vect)
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
    }
    idx_mat.close();
}

int IndexMain(int argc, char **argv)
{
    IndexWelcome();

    std::clock_t begin_time = clock();
    std::string out_dir, count_tab_path;
    size_t k_len(0);
    bool kmer_mode(false), stranded(true);

    ParseOptions(argc, argv, count_tab_path, out_dir, kmer_mode, k_len, stranded);
    PrintRunInfo(count_tab_path, out_dir, kmer_mode, k_len, stranded);

    if (!kmer_mode)
    {
        std::cerr << BOLDYELLOW << "[warning]" << RESET << " indexing in general: features are not considered as k-mers" << std::endl
                  << std::endl;
    }

    std::ifstream count_tab(count_tab_path);
    if (!count_tab.is_open())
    {
        throw std::invalid_argument("cannot open count table file: " + count_tab_path);
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = count_tab_path.find_last_of(".");
        if (pos != std::string::npos && count_tab_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(count_tab);
    std::istream kmer_count_instream(&inbuf);

    std::ofstream idx_meta(out_dir + "/idx-meta.bin"), idx_pos(out_dir + "/idx-pos.bin"), idx_mat(out_dir + "/idx-mat.bin");
    if (!idx_meta.is_open() || !idx_pos.is_open() || !idx_mat.is_open())
    {
        throw std::invalid_argument("output folder for index does not exist: " + out_dir);
    }
    ScanIndex(idx_meta, idx_pos, idx_mat, kmer_count_instream, k_len, stranded);

    idx_mat.close();
    idx_pos.close();
    idx_meta.close();
    count_tab.close();

    std::cerr << "Count table indexing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    TestIndex(out_dir + "/idx-meta.bin", out_dir + "/idx-pos.bin", out_dir + "/idx-mat.bin");

    return EXIT_SUCCESS;
}
