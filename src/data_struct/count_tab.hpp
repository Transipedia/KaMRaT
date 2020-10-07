#ifndef KAMRAT_DATASTRUCT_COUNTTAB_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>
#include <cstdint>

#include "count_tab_header.hpp"

class CountTab : public CountTabHeader
{
public:
    CountTab(const std::string &idx_file_path);

    const bool IsCountsInMem() const;
    const size_t GetTableSize() const;
    const std::string &GetIndexPath() const;

    const bool AddRowAsFields(float &row_score, const std::string &line_str, std::ofstream &idx_file);
    const bool AddRowAsString(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file);

    const float GetValue(size_t row_serial, size_t valcol_serial) const;
    const float GetStr(size_t row_serial, size_t strcol_serial) const;
    const void GetCountVect(std::vector<float> &count_vect, size_t row_serial, std::ifstream &idx_file) const;
    const float AddCountVectIfCoherent(std::vector<float> &sum_count_vect, size_t row_serial, const std::vector<float> &ref_count_vect,
                                       const std::string &dist_method, const float dist_thres, std::ifstream &idx_file) const;
    const void GetRowString(std::string &row_string, size_t row_serial, std::ifstream &idx_file, const std::string &as_scorecol) const;

    const float CalcCountDistance(size_t row_serial1, size_t row_serial2, const std::string &dist_method, std::ifstream &idx_file) const;

private:
    const std::string idx_file_path_;               // file path for count index, empty if with in-memory mode
    std::vector<std::vector<float>> value_tab_;     // value vector ------ used both in inMem or onDsk
    std::vector<std::vector<float>> count_tab_;     // count vector ------ only used in inMem
    std::vector<std::vector<std::string>> str_tab_; // string vector ----- used both in inMem or onDsk
    std::vector<size_t> index_pos_;                 // indexed position for k-mer count ------ only used in onDsk
};

using code2serial_t = std::map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_COUNTTAB_HPP