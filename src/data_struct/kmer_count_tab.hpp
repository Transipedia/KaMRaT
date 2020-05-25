#ifndef KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP
#define KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP

#include <unordered_map>
#include <vector>
#include <string>
#include <set>
#include <fstream>

#include "column_info.hpp"

class KMerCountTab
{
public:
    KMerCountTab(const std::string &mode);
    const std::string GetMode() const;
    const float AddKMerCountInMem(uint64_t kmer_code,
                                  const std::string &line_str,
                                  const ColumnInfo &column_info,
                                  const std::string &score_colname);
    const float AddKMerIndexOnDsk(uint64_t kmer_code,
                                  const std::string &line_str,
                                  const ColumnInfo &column_info,
                                  const std::string &score_colname,
                                  std::ofstream &index_file);
    const bool IsKMerExist(uint64_t kmer_code) const;
    const float GetValue(uint64_t kmer_code, size_t i_valcol) const;
    const bool GetCountInMem(std::vector<float> &count_vect, uint64_t kmer_code) const;
    const bool GetCountOnDsk(std::vector<float> &count_vect, uint64_t kmer_code, std::ifstream &index_file, size_t nb_sample) const;
    const size_t GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set) const;
    const size_t GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set, std::ifstream &index_file, size_t nb_sample) const;

private:
    const std::string mode_;                    // inMem or onDsk
    size_t nb_value_, nb_count_;                // number of value and count columns
    std::unordered_map<uint64_t, size_t> kmer_serial_;    // from k-mer unique code to line number in count table (inMem) or to position in index file (onDsk)
    std::vector<std::vector<float>> value_tab_; // k-mer value vector ------ used both in inMem or onDsk
    std::vector<std::vector<float>> count_tab_; // k-mer count vector ------ only used in inMem
    std::vector<size_t> index_pos_;             // indexed position for k-mer count ------ only used in onDsk
};

#endif //KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP