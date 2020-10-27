#ifndef KAMRAT_DATASTRUCT_FEATURETABELEM_HPP
#define KAMRAT_DATASTRUCT_FEATURETABELEM_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

class FeatureTabElem
{
public:
    FeatureTabElem(size_t nb_value, size_t nb_count, size_t nb_str,
                   std::istringstream &line_conv, std::ofstream &idx_file, const std::vector<char> &colnature_vect);

    const float GetValueAt(size_t colpos) const;
    const std::vector<float> &GetCountVect() const;
    const std::vector<std::string> &GetStrVect() const;

    const void RestoreCountVect(std::ifstream &idx_file, size_t nb_count); // for onDsk
    const void ClearCountVect();                                           // for onDsk

private:
    std::vector<float> value_vect_;     // value vector ------ used both in inMem or onDsk
    std::vector<float> count_vect_;     // count vector ------ only used in inMem
    std::vector<std::string> str_vect_; // string vector ----- used both in inMem or onDsk
    const size_t index_pos_;            // indexed position for k-mer count ------ only used in onDsk
};

using code2serial_t = std::unordered_map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number
using featuretab_t = std::vector<FeatureTabElem>;           // feature table

#endif //KAMRAT_DATASTRUCT_FEATURETABELEM_HPP