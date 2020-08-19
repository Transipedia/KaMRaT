#ifndef KAMRAT_DATASTRUCT_COUNTTAB_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>

#include "count_tab_header.hpp"

class CountTabByFields : public CountTabHeader
{
public:
    CountTabByFields(const std::string &mode);
    const std::string GetMode() const;
    const bool AddCountInMem(float &row_score, const std::string &line_str);
    const bool AddIndexOnDsk(float &row_score, const std::string &line_str, std::ofstream &index_file);
    const float GetValue(size_t row_serial, size_t valcol_serial) const;
    const float GetStr(size_t row_serial, size_t strcol_serial) const;
    const void GetCountInMem(std::vector<float> &count_vect, size_t row_serial) const;
    const void GetCountOnDsk(std::vector<float> &count_vect, size_t row_serial, std::ifstream &index_file) const;
    const size_t GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<size_t> row_serial_set) const;
    const size_t GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<size_t> row_serial_set, std::ifstream &index_file) const;

private:
    const std::string mode_;                        // inMem or onDsk
    std::vector<std::vector<float>> value_tab_;     // value vector ------ used both in inMem or onDsk
    std::vector<std::vector<float>> count_tab_;     // count vector ------ only used in inMem
    std::vector<std::vector<std::string>> str_tab_; // string vector ----- used both in inMem or onDsk
    std::vector<size_t> index_pos_;                 // indexed position for k-mer count ------ only used in onDsk
};

using code2serial_t = std::map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_COUNTTAB_HPP