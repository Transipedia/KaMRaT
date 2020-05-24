#ifndef KAMRAT_DATASTRUCT_COLUMNINFO_HPP
#define KAMRAT_DATASTRUCT_COLUMNINFO_HPP

#include <string>
#include <vector>
#include <map>

class ColumnInfo
{
public:
    ColumnInfo();
    void MakeColumnInfo(const std::string &header_line,
                        const std::string &sample_info_path);
    const std::string GetColumnName(size_t i_col) const;
    const int GetColumnLabel(size_t i_col) const;
    const char GetColumnNature(size_t i_col) const;
    const size_t GetColumnSerial(const std::string &colname) const;
    const size_t GetNbCondition() const;
    const size_t GetNbSample() const;
    const size_t GetNbColumn() const;

private:
    size_t nb_condition_, nb_sample_, nb_value_;
    std::vector<std::string> colname_vect_;        // starting from column 0
    std::map<std::string, int> colname2label_;     // non-negative for sample count (summable), -1 for value (non-summable), -2 for tag (string)
    std::map<std::string, size_t> colname2serial_; // column name to serial dictionary, serial number is assigned SEPARATELY to sample count or to value
                                                   // this attribute is useful for independent treatment of count and value in KMerCountTable object
};

#endif //KAMRAT_DATASTRUCT_COLUMNINFO_HPP