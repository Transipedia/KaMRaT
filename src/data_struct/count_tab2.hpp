#ifndef KAMRAT_DATASTRUCT_COUNTTAB2_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB2_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <unordered_map>
#include <cstdint>

#include "count_tab_header.hpp"
#include "count_tab_elem.hpp"

class CountTab2 : public CountTabHeader
{
    /// now can open/close the table file inside the class ! ///
public:
    CountTab(size_t k_len, bool stranded, const std::string &idx_file_path);

    const bool IsCountsInMem() const;
    const size_t GetTableSize() const;
    const size_t GetKLen() const;
    const bool IsStranded() const;
    const std::string &GetIndexPath() const;

    const bool AddRowAsFields(float &row_score, const std::string &line_str, std::ofstream &idx_file);
    const bool AddRowAsString(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file);

    const float GetValue(size_t row_serial, size_t valcol_serial) const;
    const float GetStr(size_t row_serial, size_t strcol_serial) const;
    const void GetCountVect(std::vector<float> &count_vect, size_t row_serial, std::ifstream &idx_file) const;
    const void EstimateMeanCountVect(std::vector<float> &avg_count_vect, const std::vector<size_t> &row_serial_vect, std::ifstream &idx_file) const;
    const void GetRowString(std::string &row_string, size_t row_serial, std::ifstream &idx_file, const std::string &as_scorecol) const;
    const float CalcCountDistance(size_t row_serial1, size_t row_serial2, const std::string &dist_method, std::ifstream &idx_file) const;

    const void ShrinkTab();

private:
    const size_t k_len_;                    // k-mer length
    const bool stranded_;                   // is k-mers stranded
    const std::string idx_file_path_;       // file path for count index, empty if with in-memory mode
    std::vector<CountTabElem> feature_tab_; // feature table
};

using code2serial_t = std::unordered_map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_COUNTTAB2_HPP