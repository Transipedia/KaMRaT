#include <vector>
#include <fstream>
#include <map>

// using code2kmer_t = std::map<uint64_t, std::pair<std::string, size_t>>;
using ftVect_t = std::vector<std::pair<std::string, size_t>>;

// uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp, called by LoadFeaturePos

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path)
{
    std::ifstream idx_meta(idx_meta_path);
    if (!idx_meta.is_open())
    {
        throw std::invalid_argument("loading index-meta failed, KaMRaT index folder not found or may be corrupted");
    }
    std::string term;
    idx_meta >> nb_smp_all >> k_len;
    if (k_len > 0)
    {
        idx_meta >> term;
        stranded = (term == "T" ? true : false);
    }
    for (size_t i(0); i <= nb_smp_all; ++i) // including the first column name
    {
        idx_meta >> term;
        colname_vect.emplace_back(term);
    }
    for (size_t i(0); i < nb_smp_all; ++i)
    {
        idx_meta >> term;
        smp_sum_vect.emplace_back(std::stod(term));
    }
    idx_meta.close();
}

void LoadPosVect(std::map<uint64_t, size_t> &code_pos_map, const std::string &idx_pos_path)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    uint64_t code;
    size_t pos;
    while (idx_pos.read(reinterpret_cast<char *>(&code), sizeof(uint64_t)) && idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        code_pos_map.insert({code, pos});
    }
    idx_pos.close();
}

// void MakeKMerMap(code2kmer_t &code_kmer_map, const std::string &idx_pos_path, const size_t k_len, const bool stranded)
// {
//     std::ifstream idx_pos(idx_pos_path);
//     if (!idx_pos.is_open())
//     {
//         throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
//     }
//     std::string term;
//     size_t pos;
//     while (idx_pos >> term >> pos)
//     {
//         size_t code = Seq2Int(term, k_len, stranded);
//         code_kmer_map.insert({code, std::make_pair(term, pos)});
//     }
//     idx_pos.close();
// }

void MakeFeatureVect(ftVect_t &feature_vect, const std::string &idx_pos_path)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    std::string term;
    size_t pos;
    while (idx_pos >> term >> pos)
    {
        feature_vect.emplace_back(std::make_pair(term, pos));
    }
    idx_pos.close();
}

const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp)
{
    count_vect.resize(nb_smp);
    idx_mat.seekg(pos);
    idx_mat.read(reinterpret_cast<char *>(&count_vect[0]), nb_smp * sizeof(float));
    return count_vect;
}

const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp,
                                       const std::vector<double> &nf_vect)
{
    GetCountVect(count_vect, idx_mat, pos, nb_smp);
    for (size_t i(0); i < nb_smp; ++i)
    {
        count_vect[i] *= nf_vect[i];
    }
    return count_vect;
}