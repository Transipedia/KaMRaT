#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <map>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "index_runinfo.hpp"
#include "index_loading.hpp"
#include "seq_coding.hpp"

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

const size_t CountColumn(const std::string &line_str)
{
    size_t nb_smp(0);
    std::istringstream conv(line_str);
    std::string term;
    // count sample number, skipping the first column
    for (conv >> term; conv >> term; ++nb_smp){}
    if (nb_smp == 0)
    {
        throw std::domain_error("input table parsing failed: sample number equals to 0");
    }
    return nb_smp;
}

void ComputeNF(std::vector<double> &nf_vect, std::istream &kmer_count_instream, const size_t nf_base)
{
    std::string line_str, term;
    std::istringstream conv;
    std::vector<float> count_vect;

    std::getline(kmer_count_instream, line_str);
    size_t nb_smp = CountColumn(line_str);
    nf_vect.resize(nb_smp, 0);
    while (std::getline(kmer_count_instream, line_str))
    {
        conv.str(line_str);
        // parse feature name and following count columns
        for (conv >> term; conv >> term; count_vect.push_back(std::stof(term))) {}
        // check if all rows have same number of columns as the header row
        if (count_vect.size() != nb_smp) 
        {
            throw std::length_error("sample numbers are not consistent: " + std::to_string(nb_smp) + " vs " + std::to_string(count_vect.size()));
        }
        for (size_t i_smp(0); i_smp < nb_smp; ++i_smp) // add count vectors together for eventual normalization
        {
            nf_vect[i_smp] += count_vect[i_smp];
        }
        conv.clear();
        count_vect.clear();
    }
    for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
    {
        nf_vect[i_smp] = nf_base / nf_vect[i_smp];
        if (nf_vect[i_smp] < 0.1)
        {
            throw std::invalid_argument("normalization factor too small (" + std::to_string(nf_vect[i_smp]) + "), please try larger base");
        }
    }
}

void IndexCount(std::ofstream &idx_pos, std::ofstream &idx_mat, const std::vector<double> &nf_vect,
                const std::string &line_str, const size_t k_len, const bool stranded, const size_t nb_smp, const bool to_norm)
{
    static std::istringstream conv(line_str);
    static std::vector<float> count_vect;
    static std::string ft_name, term;
    static std::unordered_set<uint64_t> code_set;
    static size_t ft_code(0), ft_pos;

    conv.str(line_str);
    for (conv >> ft_name; conv >> term; count_vect.push_back(std::stof(term))) // parse feature name and following count columns
    {
    }
    ft_pos = static_cast<size_t>(idx_mat.tellp());
    if (k_len > 0) // if index in k-mer mode => ft_code calculated by Seq2Int
    {
        if (k_len != ft_name.size()) // check k-mer length
        {
            throw std::length_error("feature length checking failed: length of " + ft_name + " not equal to " + std::to_string(k_len));
        }
        ft_code = Seq2Int(ft_name, k_len, stranded);
        idx_pos.write(reinterpret_cast<char *>(&ft_code), sizeof(uint64_t)); // [idx_pos] if indexing k-mer, write also k-mer code
        if (!code_set.insert(ft_code).second)
        {
            throw std::domain_error("unicity checking failed, an equivalent key already existed for feature: " + ft_name);
        }
    }
    idx_pos.write(reinterpret_cast<char *>(&ft_pos), sizeof(size_t)); // [idx_pos] feature code and feature position, ordered by code

    if (to_norm)
    {
        for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
        {
            count_vect[i_smp] *= nf_vect[i_smp];
        }
    }
    idx_mat.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(float)); // [idx_mat] feature count vector
    idx_mat << ft_name << std::endl;
    count_vect.clear();
    conv.clear();
}

void ScanIndex(std::ofstream &idx_meta, std::ofstream &idx_pos, std::ofstream &idx_mat, std::istream &kmer_count_instream,
               const std::vector<double> &nf_vect, const size_t k_len, const bool stranded, const size_t nf_base)
{
    std::string line_str;
    std::getline(kmer_count_instream, line_str); // read header row in table
    size_t nb_smp = CountColumn(line_str);

    idx_meta << nb_smp << "\t" << k_len; // [idx_meta 1] sample number, k-mer length, and strandedness if applicable
    if (k_len != 0)
    {
        idx_meta << "\t" << (stranded ? 'T' : 'F');
    }
    idx_meta << std::endl;
    idx_meta << line_str << std::endl; // [idx_meta 2] the header row
    while (std::getline(kmer_count_instream, line_str))
    {
        IndexCount(idx_pos, idx_mat, nf_vect, line_str, k_len, stranded, nb_smp, !nf_vect.empty()); // [idx_pos, idx_mat] (inside)
    }
}

void TestIndex(const std::string &idx_meta_path, const std::string &idx_pos_path, const std::string &idx_mat_path)
{
    std::ifstream idx_mat(idx_mat_path);
    std::string term;
    size_t nb_smp_all, k_len;
    bool stranded;
    std::vector<std::string> colnames_vect;
    LoadIndexMeta(nb_smp_all, k_len, stranded, colnames_vect, idx_meta_path);
    std::cout << nb_smp_all << "\t" << k_len;
    if (k_len > 0)
    {
        std::cout << "\t" << (stranded ? "T" : "F");
    }
    std::cout << std::endl;
    for (auto const &term : colnames_vect)
    {
        std::cout << term << "\t";
    }
    std::cout << std::endl;

    std::vector<size_t> pos_vect;
    LoadPosVect(pos_vect, idx_pos_path, k_len != 0);
    std::vector<float> count_vect;
    for (const size_t p : pos_vect)
    {
        GetCountVect(count_vect, idx_mat, p, nb_smp_all);
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
    std::string out_dir, count_tab_path, nf_file_path;
    size_t k_len(0), nf_base(0);
    bool stranded(true);

    ParseOptions(argc, argv, count_tab_path, out_dir, k_len, stranded, nf_base, nf_file_path);
    PrintRunInfo(count_tab_path, out_dir, k_len, stranded, nf_base, nf_file_path);
    if (0 == k_len)
    {
        std::cerr << BOLDYELLOW << "[warning]" << RESET << " indexing in general: features are not considered as k-mers" << std::endl
                  << std::endl;
    }
    if (0 == nf_base && nf_file_path.empty())
    {
        std::cerr << BOLDYELLOW << "[warning]" << RESET << " indexing without normalization" << std::endl
                  << std::endl;
    }
    else if (nf_file_path.empty())
    {
        std::cerr << BOLDYELLOW << "[warning]" << RESET << " no precomputed normalization factor given, k-mer count matrix will be scanned twice" << std::endl
                  << std::endl;
    }

    // Create and open outfiles
    std::ofstream idx_meta(out_dir + "/idx-meta.bin"), idx_pos(out_dir + "/idx-pos.bin"), idx_mat(out_dir + "/idx-mat.bin");
    if (!idx_meta.is_open() || !idx_pos.is_open() || !idx_mat.is_open())
    {
        throw std::invalid_argument("output folder for index does not exist: " + out_dir);
    }

    std::vector<double> nf_vect;
    if (nf_base > 0 && nf_file_path.empty()) // to compute NF
    {
        std::cerr << "Computing NF..." << std::endl;
        std::ifstream count_tab(count_tab_path);
        if (!count_tab.is_open())
        {
            throw std::invalid_argument("cannot open count table file: " + count_tab_path);
        }
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        if (count_tab_path.substr(count_tab_path.size() - 2) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
        inbuf.push(count_tab);
        std::istream kmer_count_instream(&inbuf);
        ComputeNF(nf_vect, kmer_count_instream, nf_base);
        count_tab.close();
    }
    else if (!nf_file_path.empty()) // to load NF
    {
        std::cerr << "Loading NF..." << std::endl;
        std::ifstream nf_file(nf_file_path);
        if (!nf_file.is_open())
        {
            throw std::invalid_argument("cannot open count NF file: " + nf_file_path);
        }
        for (double x(0); nf_file >> x; nf_vect.push_back(x))
        {
        }
        nf_file.close();
    }

    std::ifstream count_tab(count_tab_path);
    if (!count_tab.is_open())
    {
        throw std::invalid_argument("cannot open count table file: " + count_tab_path);
    }
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (count_tab_path.substr(count_tab_path.size() - 2) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(count_tab);
    std::istream kmer_count_instream(&inbuf);
    // Load and index the matrix
    ScanIndex(idx_meta, idx_pos, idx_mat, kmer_count_instream, nf_vect, k_len, stranded, nf_base);
    // Write normalization factor values to idx-meta file
    if (!nf_vect.empty())
    {
        idx_meta << nf_vect[0];
        for (size_t i(1); i < nf_vect.size(); idx_meta << "\t" << nf_vect[i])
        {
        }
    }
    
    idx_meta << std::endl;
    count_tab.close();
    idx_mat.close(), idx_pos.close(), idx_meta.close();

    // TestIndex(out_dir + "/idx-meta.bin", out_dir + "/idx-pos.bin", out_dir + "/idx-mat.bin");

    std::cerr << "Count table indexing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}
