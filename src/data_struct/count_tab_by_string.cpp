#include <sstream>
#include <iostream>

#include "count_tab_by_string.hpp"

const bool CountTabByString::IndexWithString(float &score_value, std::vector<float> &count_vect, const std::string &line_str, const size_t idx_pos)
/* ------------------------------------------------------------------------------------------- *\
    Arg:    1. string of input k-mer table line
            2. column name for score
    Value:  inserted k-mer score (a certain column's value or -1)
    Func:   1. check whether the value & count vectors to insert are coherent with the tables
            2. parse the line string
            3. insert the k-mer value & count vector to k-mer value & count table
            4. return k-mer's corresponded score value
\* ------------------------------------------------------------------------------------------- */
{
    std::vector<float> value_vect;
    std::vector<std::string> str_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, str_vect, line_str);
    input_pos_.emplace_back(idx_pos);
    if (score_pos == 0)
    {
        score_value = 0;
        return false;
    }
    else
    {
        score_value = value_vect[colserial_vect_[score_pos]];
        return true;
    }
}

const void CountTabByString::PrintFromInput(const size_t row_serial, std::ifstream &input_file, const std::string &insert_at_second_col) const
{
    size_t disk_pos = input_pos_.at(row_serial);
    std::string str_line;
    input_file.seekg(disk_pos);
    std::getline(input_file, str_line);
    std::cout << str_line.substr(0, str_line.find_first_of(" \t") + 1)
              << insert_at_second_col
              << str_line.substr(str_line.find_first_of(" \t")) << std::endl;
}
