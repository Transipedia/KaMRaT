#ifndef KAMRAT_DATASTRUCT_FEATURETAB_HPP
#define KAMRAT_DATASTRUCT_FEATURETAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <unordered_map>
#include <cstdint>

#include "feature_tab_header.hpp"
#include "feature_tab_elem.hpp"

class FeatureTab
{
public:
    FeatureTab(const std::string &idx_file_path,
               const std::string &sample_info_path,
               const std::unordered_set<std::string> &preserved_cond_tags = std::unordered_set<std::string>());

    const std::string &GetIndexPath() const;
    const size_t GetTableSize() const;

    const void MakeTable(const std::string &table_path, const std::string &repval_name);
    // const void MakeHeaderRow(const std::string &header_str, const std::string &repval_name);
    // const void AddRowAsFields(float &row_score, const std::string &line_str, std::ofstream &idx_file);
    // const void IndexStringRow(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file);

    const float GetValue(size_t row_serial, size_t valcol_serial) const;
    const float GetRepValue(size_t row_serial) const;
    const float GetStr(size_t row_serial, size_t strcol_serial) const;
    const void GetCountVect(std::vector<float> &count_vect, size_t row_serial, std::ifstream &idx_file);
    const void EstimateMeanCountVect(std::vector<float> &avg_count_vect, const std::vector<size_t> &row_serial_vect, std::ifstream &idx_file) const;
    const void GetRowString(std::string &row_string, size_t row_serial, std::ifstream &idx_file, const std::string &as_scorecol) const;
    const float CalcCountDistance(size_t row_serial1, size_t row_serial2, const std::string &dist_method, std::ifstream &idx_file) const;

    const void ShrinkTab();

private:
    size_t nb_row_, repval_colpos_;           // different column number, rep-value position, and table size
    const std::string idx_file_path_;         // file path for count index, empty if with in-memory mode
    FeatureTabHeader feature_tab_header_;     // header line info
    std::vector<FeatureTabElem> feature_tab_; // feature table
};

using code2serial_t = std::unordered_map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_FEATURETAB_HPP