#ifndef KAMRAT_DATASTRUCT_COUNTTABHEADER_HPP
#define KAMRAT_DATASTRUCT_COUNTTABHEADER_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>

class CountTabHeader
{
public:
    CountTabHeader();
    const std::string GetMode() const;
    //----- header info -----//
    const void MakeColumnInfo(const std::string &header_line, const std::string &sample_info_path, const std::string &score_colname);
    const std::string GetColName(size_t i_col) const;
    const char GetColNature(size_t i_col) const;
    const size_t GetColSerial(size_t i_col) const;
    const size_t GetSmpLabel(size_t i_col) const;
    const std::vector<size_t> &GetSmpLabels() const;
    const size_t GetNbCondition() const;
    const size_t GetNbValue() const;
    const size_t GetNbCount() const;
    const size_t GetNbColumn() const;

protected:
    //----- header info -----//
    size_t nb_cond_, nb_value_, nb_count_, nb_str_; // number of conditions and value/count/string columns
    std::vector<std::string> colname_vect_;         // column name vector parsed from the header line
    std::vector<char> colnature_vect_;              // column nature vector, f=feature s=sample v=value +=score_value
    std::vector<size_t> colserial_vect_;            // column serial vector, serial number is assigned SEPARATELY to count or to value
    std::vector<size_t> smplabel_vect_;             // sample label vector parsed from the header line and sample-info file
};

#endif //KAMRAT_DATASTRUCT_COUNTTABHEADER_HPP