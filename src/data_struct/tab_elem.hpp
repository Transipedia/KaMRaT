#ifndef KAMRAT_DATASTRUCT_TABELEM_HPP
#define KAMRAT_DATASTRUCT_TABELEM_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "tab_header.hpp"

class TabElem
{
public:
    TabElem(std::istringstream &line_conv, std::ofstream &idx_file, const TabHeader &tab_header);                                 // Index the table row
    TabElem(std::istringstream &line_conv, std::vector<float> &count_vect, std::ofstream &idx_file, const TabHeader &tab_header); // Index with given vector, string, and value

    /* Reorganize header string: firstly non-count columns, then count columns */
    const std::string &MakeOutputRowStr(std::string &row_str, std::ifstream &idx_file, size_t nb_count);                      // with indexed info
    const std::string &MakeOutputRowStr(std::string &row_str, const std::vector<float> &count_vect, std::ifstream &idx_file); // with given count vector

    const float GetValue() const;                                                                                           // Get row value
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, size_t nb_count) const; // Get count vector
    const std::string &GetNonCountStr(std::string &non_count_str, std::ifstream &idx_file, size_t nb_count) const;          // Get non-count string
    const size_t GetIdxPos() const;                                                                                         // Get index position

protected:
    float value_;                   // representative value or score
    std::vector<float> count_vect_; // count vector ------ only used in inMem
    std::string non_count_str_;     // non-count string ------ only used in inMem
    size_t idx_pos_;                // row index position ------ only used in onDsk
};

using countTab_t = std::vector<TabElem>;                    // feature table
using code2serial_t = std::unordered_map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_TABELEM_HPP