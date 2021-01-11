#include <iostream>
#include <cmath>

#include "kmer_elem.hpp"

const float TransformValue(const float raw_value, const RepModeCode rep_mode_code)
{
    static size_t input_serial = 0;
    input_serial++;
    switch (rep_mode_code)
    {
    case RepModeCode::kMin:
        return raw_value;
    case RepModeCode::kMinAbs:
        return fabs(raw_value);
    case RepModeCode::kMax:
        return (-raw_value);
    case RepModeCode::kMaxAbs:
        return -fabs(raw_value);
    default: // RepModeCode::kInputOrder
        return input_serial;
    }
}

KMerElem::KMerElem(const float value, std::vector<float> &count_vect, const std::string &value_str,
                   std::ofstream &idx_file, const RepModeCode rep_mode_code)
    : idx_pos_(static_cast<size_t>(idx_file.tellp())),
      rep_value_(TransformValue(value, rep_mode_code))
{
    idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(float)); // first index counts columns
    idx_file << value_str << std::endl;                                                          // then index non-count columns
}

const float KMerElem::GetRepValue() const noexcept
{
    return rep_value_;
}

const std::vector<float> &KMerElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count) const
{
    count_vect.resize(nb_count);
    idx_file.seekg(idx_pos_);
    idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(float));
    return count_vect;
}

const float KMerElem::GetCountAt(std::ifstream &idx_file, const size_t i_smp) const
{
    static float count;
    idx_file.seekg(idx_pos_ + i_smp * sizeof(float));
    idx_file.read(reinterpret_cast<char *>(&count), sizeof(float));
    return count;
}

const std::string &&KMerElem::GetValueStr(std::ifstream &idx_file, const size_t nb_count) const
{
    static std::string value_str;
    value_str.clear();
    idx_file.seekg(idx_pos_ + nb_count * sizeof(float));
    std::getline(idx_file, value_str);
    return std::move(value_str);
}