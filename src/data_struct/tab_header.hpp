#ifndef KAMRAT_DATASTRUCT_TABHEADER_HPP
#define KAMRAT_DATASTRUCT_TABHEADER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

/* ============= Notation for column nature ============= *\
 *     0                indicates "not a sample"          *
 *     positive value   corresponds to condition label    *
\* ====================================================== */

class TabHeader
{
public:
    TabHeader();                                 // all columns except the first are about sample count
    TabHeader(const std::string &smp_info_path); // to indicate which columns are about sample count

    const void MakeColumns(std::istringstream &line_conv, const std::string &rep_colname); // Make columns from header row string

    const size_t GetNbCondition() const; // Number of conditions
    const size_t GetNbCol() const;       // Number of columns
    const size_t GetNbCount() const;     // Number of count columns
    const size_t GetRepColPos() const;   // Position of representative value column

    const std::string &GetColNameAt(size_t i_col) const; // Get column name at given column number
    const size_t GetColNatAt(size_t i_col) const;        // Get column nature at given column number
    const bool IsColCount(size_t i_col) const;           // Is given column number a count column

    const std::vector<size_t> &GetCondiSmpNbVect() const;          // Get number of samples in each condition

    const double ParseRowStr(std::vector<float> &count_vect,
                             std::string &non_count_str,
                             std::istringstream &line_conv) const; // Parse the table row string

private:
    size_t rep_colpos_;                                              // position of column indicating representative or score value
    std::vector<std::string> condi_name_vect_;                       // vector of condition names
    std::unordered_map<std::string, size_t> smp2lab_;                // sample name to label, condition name to label
    std::vector<size_t> colnat_nb_vect_;                             // vector of sample number in each type of column
    std::vector<std::pair<std::string, uint64_t>> col_namenat_vect_; // vector of column name-nature pairs
};

#endif //KAMRAT_DATASTRUCT_TABHEADER_HPP