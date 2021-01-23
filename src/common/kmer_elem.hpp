#ifndef KAMRAT_COMMON_KMERELEM_HPP
#define KAMRAT_COMMON_KMERELEM_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <memory>

#include "tab_header.hpp"

/* ======================================= *\
 * Used by kamratMerge and contigEvaluator *
\* ======================================= */

enum RepModeCode
{
    kMin = 0,
    kMinAbs,
    kMax,
    kMaxAbs,
    kInputOrder
};

class KMerElem
{
public:
    KMerElem(float value, std::vector<float> &count_vect, const std::string &value_str,
             std::ofstream &idx_file, const RepModeCode rep_mode_code = RepModeCode::kInputOrder); // Index with given elements

    const float GetRepValue() const noexcept;                                                                               // Get row value
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, size_t nb_count) const; // Get count vector
    const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, size_t nb_count,
                                           const std::vector<double> &nf_vect) const; // Get count vector and do normalization
    const float GetCountAt(std::ifstream &idx_file, size_t i_smp) const;              // Get one count
    const std::string &&GetValueStr(std::ifstream &idx_file, size_t nb_count) const;  // Get value string

protected:
    const size_t idx_pos_;  // row index position ------ only used in onDsk
    const float rep_value_; // representative value
};

using kMerTab_t = std::vector<KMerElem>;                                     // feature table
using code2kmer_t = std::unordered_map<uint64_t, KMerElem>; // external dictionary to link k-mer with k-mer element

#endif //KAMRAT_COMMON_KMERELEM_HPP