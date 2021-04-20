#ifndef KAMRAT_DATASTRUCT_FEATUREELEM_HPP
#define KAMRAT_DATASTRUCT_FEATUREELEM_HPP

#include <vector>
#include <fstream>
#include <map>
#include <memory>

class FeatureElem
{
public:
    FeatureElem(const std::string &ft_name, size_t ft_pos);
    const std::string &GetName() const;
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, size_t nb_smp) const;
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                           std::ifstream &idx_mat, size_t pos, size_t nb_smp,
                                           const std::vector<double> &nf_vect) const;

private:
    const std::string ft_name_;
    const size_t ft_pos_;
};

using code2kmer_t = std::map<uint64_t, std::unique_ptr<FeatureElem>>;
using ftVect_t = std::vector<std::unique_ptr<FeatureElem>>;

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path);
void MakeKMerMap(code2kmer_t &code_kmer_map, const std::string &idx_pos_path, size_t k_len, bool stranded);
void MakeFeatureMap(ftVect_t &ft_vect, const std::string &idx_pos_path);

#endif //KAMRAT_DATASTRUCT_FEATUREELEM_HPP