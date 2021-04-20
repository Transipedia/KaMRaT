#ifndef KAMRAT_COMMON_KMerCOUNTTAB_HPP
#define KAMRAT_COMMON_KMerCOUNTTAB_HPP

#include <vector>
#include <string>
#include <fstream>
#include <map>

class KMerCountTab
{
public:
    KMerCountTab(std::vector<std::pair<std::string, size_t>> &kmer_vect, const std::string &idx_info_path);

    const size_t GetSampleNumber() const;
    const size_t GetKLen() const;
    const bool IsStranded() const;
    const std::map<uint64_t, size_t> &GetKMerSerialMap() const;
    // const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, size_t serial) const;
    // const std::vector<double> &GetSumVect(std::vector<double> &sum_vect, std::ifstream &idx_mat) const;
    // const float GetCountAt(std::ifstream &idx_file, size_t i_smp) const;              // Get one count

private:
    size_t nb_smp_, k_len_;                      // sample number, k-mer length
    bool stranded_;                              // k-mer strandedness if applicable
    std::vector<std::string> colname_vect_;      // column names: nb_smp_ + 1
    std::vector<size_t> pos_vect_;               // k-mer position in indexed matrix
    std::map<uint64_t, size_t> kmer_serial_map_; // k-mer serial
};

#endif //KAMRAT_COMMON_KMerCOUNTTAB_HPP