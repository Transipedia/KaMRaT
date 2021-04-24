#include <vector>
#include <fstream>
#include <map>
#include <unordered_map>
#include <algorithm>

// using code2kmer_t = std::map<uint64_t, std::pair<std::string, size_t>>;
using ftVect_t = std::vector<std::pair<std::string, size_t>>;

uint64_t Seq2Int(const std::string &seq, size_t k_length, bool stranded); // in utils/seq_coding.cpp
uint64_t GetRC(const uint64_t code, size_t k_length);                     // in utils/seq_coding.cpp
uint64_t NextCode(uint64_t code, size_t k_length, char new_nuc);          // in utils/seq_coding.cpp

const double CalcVectMedian(const std::vector<float> &x)
{
    static std::vector<float> x_tmp;
    x_tmp = x;
    size_t n_elem = x_tmp.size();
    auto mid_elem = x_tmp.begin() + n_elem / 2;
    std::nth_element(x_tmp.begin(), mid_elem, x_tmp.end());
    if (n_elem % 2 != 0)
    {
        return *mid_elem;
    }
    else
    {
        return 0.5 * (*mid_elem + *(std::max_element(x_tmp.begin(), mid_elem)));
    }
}

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

void LoadPosVect(std::vector<size_t> &pos_vect, const std::string &idx_pos_path)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    size_t pos;
    while (idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        pos_vect.emplace_back(pos);
    }
    idx_pos.close();
}

void LoadCodePosMap(std::map<uint64_t, size_t> &code_pos_map, const std::string &idx_pos_path)
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

void LoadCodePosValMap(std::map<uint64_t, std::pair<size_t, float>> &code_posval_map, std::unordered_map<uint64_t, float> &sel_code_val_map,
                       const std::string &idx_pos_path)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    uint64_t code;
    size_t pos;
    std::unordered_map<uint64_t, float>::const_iterator it;
    const bool empty_sel_list = sel_code_val_map.empty();
    while (idx_pos.read(reinterpret_cast<char *>(&code), sizeof(uint64_t)) && idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        it = sel_code_val_map.find(code);
        if (it != sel_code_val_map.cend()) // only load code in the given list
        {
            code_posval_map.insert({code, std::make_pair(pos, it->second)});
            sel_code_val_map.erase(it);
        }
        else if (empty_sel_list) // if the given list is empty, then load all k-mers in index
        {
            code_posval_map.insert({code, std::make_pair(pos, 0)});
        }
    }
    if (!sel_code_val_map.empty())
    {
        throw std::domain_error(std::to_string(sel_code_val_map.size()) + " k-mer(s) in the selection file do not exist in index");
    }
    idx_pos.close();
}

const std::string &GetTagSeq(std::string &tag_str, std::ifstream &idx_mat, const size_t pos, const size_t nb_smp)
{
    idx_mat.seekg(pos + nb_smp * sizeof(float));
    idx_mat >> tag_str;
    return tag_str;
}

const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp)
{
    count_vect.resize(nb_smp);
    idx_mat.seekg(pos);
    idx_mat.read(reinterpret_cast<char *>(&count_vect[0]), nb_smp * sizeof(float));
    return count_vect;
}

const std::vector<float> &GetMeanCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                           const std::string &seq, const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                                           const size_t k_len, const bool stranded)
{
    static std::vector<float> count_vect_x;

    size_t seq_len = seq.size(), nb_mem_kmer(0);
    uint64_t kmer_code = Seq2Int(seq.substr(0, k_len), k_len, true);
    std::map<uint64_t, std::pair<size_t, float>>::const_iterator it;
    count_vect.assign(nb_smp, 0);
    for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
    {
        it = code_posval_map.find(kmer_code);
        if (it == code_posval_map.cend() && !stranded)
        {
            it = code_posval_map.find(GetRC(kmer_code, k_len));
        }
        if (it != code_posval_map.cend())
        {
            GetCountVect(count_vect_x, idx_mat, it->second.first, nb_smp);
            for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
            {
                count_vect[i_smp] += count_vect_x[i_smp];
            }
            nb_mem_kmer++;
        }
        if (start_pos + k_len < seq_len)
        {
            kmer_code = NextCode(kmer_code, k_len, seq[start_pos + k_len]);
        }
    }
    for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
    {
        count_vect[i_smp] /= nb_mem_kmer;
    }

    return count_vect;
}

const std::vector<float> &GetMedianCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                             const std::string &seq, const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                                             const size_t k_len, const bool stranded)
{
    static std::vector<std::vector<float>> mem_kmer_counts(nb_smp);
    static std::vector<float> count_vect_x;

    size_t seq_len = seq.size();
    uint64_t kmer_code = Seq2Int(seq.substr(0, k_len), k_len, true);
    std::map<uint64_t, std::pair<size_t, float>>::const_iterator it;
    for (size_t start_pos(0); start_pos < seq_len - k_len + 1; ++start_pos)
    {
        it = code_posval_map.find(kmer_code);
        if (it == code_posval_map.cend() && !stranded)
        {
            it = code_posval_map.find(GetRC(kmer_code, k_len));
        }
        if (it != code_posval_map.cend())
        {
            GetCountVect(count_vect_x, idx_mat, it->second.first, nb_smp);
            for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
            {
                mem_kmer_counts[i_smp].emplace_back(count_vect_x[i_smp]);
            }
        }
        if (start_pos + k_len < seq_len)
        {
            kmer_code = NextCode(kmer_code, k_len, seq[start_pos + k_len]);
        }
    }
    for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
    {
        count_vect.push_back(CalcVectMedian(mem_kmer_counts[i_smp]));
    }

    return count_vect;
}