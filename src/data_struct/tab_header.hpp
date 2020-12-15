#ifndef KAMRAT_DATASTRUCT_TABHEADER_HPP
#define KAMRAT_DATASTRUCT_TABHEADER_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

/* ====== Notation for column nature ====== *\
 *     A-Z:    sample labels                *
 *     v:      values                       *
 *     s:      strings                      *
\* ======================================== */

class TabHeader
{
public:
    TabHeader();                                 // all columns except the first are about sample count
    TabHeader(const std::string &smp_info_path); // to indicate which columns are about sample count

    const void MakeColumnInfo(std::istringstream &line_conv, const std::string &rep_colname); // Parse the header row string

    const std::string &MakeOutputHeaderStr(std::string &header) const; // Reorganize header string: firstly non-count columns, then count columns

    const size_t GetNbCondition() const;                            // Number of conditions
    const size_t GetConditionLabel(const std::string &condi) const; // Get condition label from condition name
    const size_t GetNbCol() const;                                  // Number of columns
    const size_t GetNbCount() const;                                // Number of count columns
    const size_t GetRepColPos() const;                              // Position of representative value column
    const std::string &GetColNameAt(size_t i) const;                // Get column name at given column number
    const bool IsCount(size_t i_col) const;                         // Is given column number a count column

    const double ParseRowStr(std::vector<float> &count_vect, std::string &non_count_str,
                             std::istringstream &line_conv) const; // Parse the table row string
    const std::vector<size_t> &GetSmpLabels() const;               // Get sample labels

    const void PrintSmp2Lab() const; // Print sample-label map, for debug

protected:
    size_t nb_count_, rep_colpos_, nb_condi_;                     // number of values, counts, strings, columns, conditions
    std::unordered_map<std::string, size_t> smp2lab_, condi2lab_; // sample name to label, condition name to label
    std::vector<std::string> colname_vect_;                       // column name vector parsed from the header line
    std::vector<bool> is_count_;                                  // 0 means not a sample, positive number indicates condition
    std::vector<size_t> sample_labels_;                           // labels of sample columns
};

#endif //KAMRAT_DATASTRUCT_TABHEADER_HPP