#ifndef KAMRAT_UTILS_COLUMNINFO_HPP
#define KAMRAT_UTILS_COLUMNINFO_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

class ColumnInfo
{
public:
    ColumnInfo();
    void MakeColumnInfo(const std::string &header_line,
                        const std::string &sample_info_path,
                        const std::string &rep_col_name);
    const char GetColumnNature(size_t i_col) const;
    const unsigned int GetNbCondition() const;
    const unsigned int GetNbSample() const;
    void PrintSampleNames(std::ostream &out_s, const std::string &leading_string) const;

private:
    unsigned int nb_condition_;
    unsigned int nb_sample_;
    std::vector<std::string> col_name_vect_;
    std::vector<int> col_nat_vect_; // starting from column 0
};

#endif //KAMRAT_UTILS_COLUMNINFO_HPP