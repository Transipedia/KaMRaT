#ifndef KAMRAT_DATASTRUCT_FEATURETABELEM_HPP
#define KAMRAT_DATASTRUCT_FEATURETABELEM_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <unordered_map>
#include <cstdint>

class FeatureTabElem
{
public:
    FeatureTabElem(size_t nb_value, size_t nb_count, size_t nb_str,
                   std::istringstream &line_conv, bool count_on_disk, const std::vector<char> &colnature_vect);

    const std::vector<float> &GetValueVect() const;
    const std::vector<float> &GetCountVect() const;
    const std::vector<std::string> &GetStrVect() const;
    const size_t GetIndexPos() const;

    const void RestoreCountVect(std::ifstream &idx_file, size_t nb_count); // for onDsk

private:
    std::vector<float> value_vect_;     // value vector ------ used both in inMem or onDsk
    std::vector<float> count_vect_;     // count vector ------ only used in inMem
    std::vector<std::string> str_vect_; // string vector ----- used both in inMem or onDsk
    size_t index_pos_;                  // indexed position for k-mer count ------ only used in onDsk
};

#endif //KAMRAT_DATASTRUCT_FEATURETABELEM_HPP