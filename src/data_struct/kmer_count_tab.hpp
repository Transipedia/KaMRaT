#ifndef KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP
#define KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP

#include <map>
#include <vector>
#include <set>
#include <string>
#include <fstream>

#include "column_info.hpp"

class KMerCountTab
{
public:
    KMerCountTab(const std::string &mode, const ColumnInfo &col_info);
    const bool AddKMerCountInMem(uint64_t kmer_code, const std::string &line_str);                            // return whether adding is succeeded
    const bool AddKMerIndexOnDsk(uint64_t kmer_code, const std::string &line_str, std::ofstream &index_file); // return whether adding is succeeded
    const bool IsKMerExist(uint64_t kmer_code) const;
    const bool GetValue(float &value, uint64_t kmer_code, const std::string &value_name) const;
    const bool GetCountInMem(std::vector<float> &count_vect, uint64_t kmer_code) const;
    const bool GetCountOnDsk(std::vector<float> &count_vect, uint64_t kmer_code, std::ifstream &index_file) const;
    const bool GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set) const;
    const bool GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<uint64_t> kmer_code_set, std::ifstream &index_file) const;

private:
    const std::string mode_;                        // inMem or onDsk
    const ColumnInfo col_info_;                     // names of info columns
    std::map<std::string, size_t> value_name_dict_; // from value's column name to column number in value table
    std::map<uint64_t, size_t> kmer_serial_;        // from k-mer unique code to line number in count table (inMem) or to position in index file (onDsk)
    std::vector<std::vector<float>> value_tab_;     // k-mer value vector (used both in inMem or onDsk)
    std::vector<std::vector<float>> count_tab_;     // k-mer count vector (only used in inMem)
    std::vector<size_t> index_pos_;                 // indexed position for k-mer count (only used in onDsk)
};

#endif //KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP