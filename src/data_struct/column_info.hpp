#ifndef KAMRAT_UTILS_COLUMNINFO_HPP
#define KAMRAT_UTILS_COLUMNINFO_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

class ColumnInfo
{
public:
    ColumnInfo(const std::string &header_line,
               const std::string &sample_info_path,
               const std::string &rep_col_name);
    const std::string GetColumnName(size_t i_col) const;
    const char GetColumnNature(size_t i_col) const;
    const unsigned int GetNbCondition() const;
    const unsigned int GetNbSample() const;
    const unsigned int GetNbColumn() const;

private:
    unsigned int nb_condition_;
    unsigned int nb_sample_;
    std::vector<std::string> col_name_vect_; // starting from column 0
    std::vector<int> col_nat_vect_;          // starting from column 0
};

#endif //KAMRAT_UTILS_COLUMNINFO_HPP