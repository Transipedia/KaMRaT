#ifndef KAMRAT_DATASTRUCT_COUNTTAB_HPP
#define KAMRAT_DATASTRUCT_COUNTTAB_HPP

#include <vector>
#include <string>
#include <fstream>

#include "count_tab_header.hpp"

class CountTabByString : public CountTabHeader
{
public:
    const bool IndexWithString(float &score_value, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &index_file);
    const void PrintFromInput(size_t row_serial, std::istream &index_file, const std::string &insert_at_second_col) const;

private:
    std::vector<size_t> index_pos_;    // indexed position for k-mer count ------ only used in onDsk
};

#endif //KAMRAT_DATASTRUCT_COUNTTAB_HPP