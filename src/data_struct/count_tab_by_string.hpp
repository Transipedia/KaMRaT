#ifndef KAMRAT_DATASTRUCT_COUNTTAB_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>

#include "count_tab_header.hpp"

class CountTabByString : public CountTabHeader
{
public:
    const bool IndexWithString(float &score_value, std::vector<float> &count_vect, const std::string &line_str, size_t idx_pos);
    const void PrintFromInput(size_t row_serial, std::ifstream &index_file, const std::string &insert_at_second_col) const;

private:
    std::vector<size_t> input_pos_;             // indexed position for k-mer count ------ only used in onDsk
};

#endif //KAMRAT_DATASTRUCT_COUNTTAB_HPP