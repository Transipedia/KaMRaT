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
    TabElem(float value, std::vector<float> &count_vect, const std::string &value_str, std::ofstream &idx_file); // Index with given elements

    const float GetRepValue() const;                                                                                        // Get row value
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, size_t nb_count) const; // Get count vector
    const std::string &GetValueStr(std::string &non_count_str, std::ifstream &idx_file, size_t nb_count) const;             // Get value string

protected:
    const size_t idx_pos_;  // row index position ------ only used in onDsk
    const float rep_value_; // representative value
};

using countTab_t = std::vector<TabElem>;                    // feature table
using code2serial_t = std::unordered_map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_TABELEM_HPP