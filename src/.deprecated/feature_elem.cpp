#include "feature_elem.hpp"

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp, called by LoadFeaturePos

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

void MakeKMerMap(code2kmer_t &code_kmer_map, const std::string &idx_pos_path, const size_t k_len, const bool stranded)
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
        size_t code = Seq2Int(term, k_len, stranded);
        code_kmer_map.insert({code, std::make_unique<FeatureElem>(term, pos)});
    }
    idx_pos.close();
}

void MakeFeatureMap(ftVect_t &feature_vect, const std::string &idx_pos_path)
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
        feature_vect.emplace_back(std::make_unique<FeatureElem>(term, pos));
    }
    idx_pos.close();
}

FeatureElem::FeatureElem(const std::string &ft_name, const size_t ft_pos)
    : ft_name_(ft_name), ft_pos_(ft_pos)
{
}

const std::string &FeatureElem::GetName() const
{
    return ft_name_;
}

const std::vector<float> &FeatureElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp) const
{
    count_vect.resize(nb_smp);
    idx_mat.seekg(ft_pos_);
    idx_mat.read(reinterpret_cast<char *>(&count_vect[0]), nb_smp * sizeof(float));
    return count_vect;
}

const std::vector<float> &FeatureElem::GetCountVect(std::vector<float> &count_vect,
                                                    std::ifstream &idx_mat, const size_t pos, const size_t nb_smp,
                                                    const std::vector<double> &nf_vect) const
{
    GetCountVect(count_vect, idx_mat, nb_smp);
    for (size_t i(0); i < nb_smp; ++i)
    {
        count_vect[i] *= nf_vect[i];
    }
    return count_vect;
}