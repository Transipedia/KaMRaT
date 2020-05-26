#ifndef KAMRAT_DATASTRUCT_COLUMNINFO_HPP
#define KAMRAT_DATASTRUCT_COLUMNINFO_HPP

#include <string>
#include <vector>

class ColumnInfo
{
public:
    ColumnInfo();
    void MakeColumnInfo(const std::string &header_line,
                        const std::string &sample_info_path,
                        const std::string &score_colname);
    const std::string GetColumnName(size_t i_col) const;
    const int GetColLabel(size_t i_col) const;
    const char GetColNature(size_t i_col) const;
    const size_t GetColSerial(size_t i_col) const;
    const size_t GetNbCondition() const;
    const size_t GetNbValue() const;
    const size_t GetNbCount() const;
    const size_t GetNbColumn() const;

private:
    size_t nb_condition_, nb_value_, nb_count_; // number of condition, value, and sample count
    std::vector<std::string> colname_vect_;     // starting from column 0
    std::vector<int> collabel_vect_;            // non-negative for sample count (summable), -1 for value (non-summable), -2 for special value, -3 for tag
    std::vector<size_t> colserial_vect_;        // column name to serial dictionary, serial number is assigned SEPARATELY to sample count or to value
                                                // this attribute is useful for independent treatment of count and value in KMerCountTable object
};

#endif //KAMRAT_DATASTRUCT_COLUMNINFO_HPP