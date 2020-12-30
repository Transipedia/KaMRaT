#include <iostream>

#include "tab_elem.hpp"

TabElem::TabElem(const float value, std::vector<float> &count_vect, const std::string &value_str, std::ofstream &idx_file)
    : idx_pos_(static_cast<size_t>(idx_file.tellp())),
      rep_value_(value)
{
    idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(count_vect[0])); // first index counts columns
    idx_file << value_str << std::endl;                                                                  // then index non-count columns
}

const float TabElem::GetRepValue() const noexcept
{
    return rep_value_;
}

const std::vector<float> &TabElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count) const
{
    count_vect.resize(nb_count);
    idx_file.seekg(idx_pos_);
    idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(float));
    return count_vect;
}

const float TabElem::GetCountAt(std::ifstream &idx_file, const size_t i_smp) const
{
    static float count;
    idx_file.seekg(idx_pos_ + i_smp * sizeof(float));
    idx_file.read(reinterpret_cast<char *>(&count), sizeof(float));
    return count;
}

const std::string &&TabElem::GetValueStr(std::ifstream &idx_file, const size_t nb_count) const
{
    static std::string value_str;
    value_str.clear();
    idx_file.seekg(idx_pos_ + nb_count * sizeof(float));
    std::getline(idx_file, value_str);
    return std::move(value_str);
}