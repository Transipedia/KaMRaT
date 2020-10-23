#ifndef KAMRAT_DATASTRUCT_FEATURETABHEADER_HPP
#define KAMRAT_DATASTRUCT_FEATURETABHEADER_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

/* ====== Notation for column nature ====== *\
 *     A-Z:    sample labels                *
 *     v:      values                       *
 *     s:      strings                      *
\* ======================================== */

class FeatureTabHeader
{
public:
    FeatureTabHeader(const std::string &sample_info_path, const std::unordered_set<std::string> &preserved_cond_tags);

    const size_t MakeColumnInfo(const std::string &header_line, const std::string &rep_colname);

    const size_t GetNbValue() const;
    const size_t GetNbCount() const;
    const size_t GetNbStr() const;
    const size_t GetNbCol() const;
    const size_t GetNbCondition() const;

    const std::vector<std::string> &GetColNameVect() const;
    const std::vector<char> &GetColNatureVect() const;
    const std::vector<size_t> &GetColSerialVect() const;

    const bool IsSample(size_t i_col) const;
    const void ParseSmpLabels(std::vector<size_t> &smp_labels);

protected:
    size_t nb_value_, nb_count_, nb_str_, nb_col_, nb_condi_;   // number of values, counts, strings, columns, conditions
    std::unordered_map<std::string, char> smp2lab_, condi2lab_; // sample name to label, condition name to label

    std::vector<std::string> colname_vect_; // column name vector parsed from the header line
    std::vector<char> colnature_vect_;      // column nature vector, cf. above for nature code
    std::vector<size_t> colserial_vect_;    // column serial vector, serial number is assigned SEPARATELY to count or to value
};

#endif //KAMRAT_DATASTRUCT_FEATURETABHEADER_HPP