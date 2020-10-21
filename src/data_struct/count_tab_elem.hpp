#ifndef KAMRAT_DATASTRUCT_COUNTTABELEM_HPP
#define KAMRAT_DATASTRUCT_COUNTTABELEM_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <unordered_map>
#include <cstdint>

#include "count_tab_header.hpp"

class CountTabElem
{
public:
    CountTabElem(size_t nb_value, size_t nb_count, size_t nb_str);

    const bool MakeRowByFields(float &row_score, const std::string &line_str);                                                          // inMem if empty index file
    const bool MakeRowByString(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file); // always onDsk

    const float GetValueAt(size_t valcol_serial) const;
    const float GetStrAt(size_t strcol_serial) const;
    const std::vector<float> &GetCountVect(bool get_by_move) const;
    const size_t GetIndexPos() const;

    const void RestoreCountVect(std::ifstream &idx_file); // for onDsk

private:
    std::vector<float> value_vect_;     // value vector ------ used both in inMem or onDsk
    std::vector<float> count_vect_;     // count vector ------ only used in inMem
    std::vector<std::string> str_vect_; // string vector ----- used both in inMem or onDsk
    size_t index_pos_;                  // indexed position for k-mer count ------ only used in onDsk
};

#endif //KAMRAT_DATASTRUCT_COUNTTABELEM_HPP