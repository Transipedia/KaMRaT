#ifndef KAMRAT_DATASTRUCT_COUNTTAB_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>

#include "count_tab_header.hpp"

class CountTab : public CountTabHeader
{
public:
    CountTab(const std::string &mode);
    const std::string GetMode() const;
    const float AddCountInMem(const std::string &line_str);
    const float AddIndexOnDsk(const std::string &line_str, std::ofstream &index_file);
    const bool IndexWithString(float &score_value, std::vector<float> &count_vect, const std::string &line_str, size_t idx_pos);
    const float GetValue(size_t kmer_serial, size_t valcol_serial) const;
    const void GetCountInMem(std::vector<float> &count_vect, size_t kmer_serial) const;
    const void GetCountOnDsk(std::vector<float> &count_vect, size_t kmer_serial, std::ifstream &index_file) const;
    const size_t GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set) const;
    const size_t GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set, std::ifstream &index_file) const;
    const void PrintFromIndex(size_t row_serial, std::ifstream &index_file, bool with_first_col) const;

private:
    const std::string mode_; // inMem or onDsk
    std::vector<std::vector<float>> value_tab_; // k-mer value vector ------ used both in inMem or onDsk
    std::vector<std::vector<float>> count_tab_; // k-mer count vector ------ only used in inMem
    std::vector<size_t> index_pos_;             // indexed position for k-mer count ------ only used in onDsk
};

using code2serial_t = std::map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_COUNTTAB_HPP