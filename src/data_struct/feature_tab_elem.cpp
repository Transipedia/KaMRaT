#include <iostream>

#include "feature_tab_elem.hpp"
#include "statistics.hpp"

FeatureTabElem::FeatureTabElem(const size_t nb_value, const size_t nb_count, const size_t nb_str,
                               std::istringstream &line_conv, std::ofstream &idx_file, const std::vector<char> &colnature_vect)
    : index_pos_(idx_file.is_open() ? idx_file.tellp() : 0)
{
    value_vect_.reserve(nb_value);
    count_vect_.reserve(nb_count);
    str_vect_.reserve(nb_str);

    static std::string term;
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    for (unsigned int i(1); i < colnature_vect.size(); ++i)
    {
        if (colnature_vect[i] == 'v') // v for values
        {
            value_vect_.emplace_back(std::move(std::stof(term)));
        }
        else if (colnature_vect[i] >= 'A' && colnature_vect[i] <= 'Z') // A-Z for sample conditions, maximally support 26 conditions
        {
            count_vect_.emplace_back(std::move(std::stof(term)));
        }
        else if (colnature_vect[i] == 's') // s for strings
        {
            str_vect_.emplace_back(std::move(term));
        }
        else
        {
            throw std::invalid_argument("unknown column nature code: " + std::to_string(colnature_vect[i]));
        }
    }
    if (line_conv >> term) // should not happen, for debug
    {
        throw std::domain_error("parsing string line to fields failed");
    }
    if (idx_file.is_open())
    {
        idx_file.write(reinterpret_cast<char *>(&count_vect_[0]), count_vect_.size() * sizeof(count_vect_[0]));
        count_vect_.clear();
    }
}

const float FeatureTabElem::GetValueAt(const size_t colpos) const
{
    return value_vect_[colpos];
}

const std::vector<float> &FeatureTabElem::GetCountVect() const
{
    if (count_vect_.empty()) // should not happen, for debug
    {
        throw std::domain_error("query empty count vector");
    }
    return count_vect_;
}

const std::vector<std::string> &FeatureTabElem::GetStrVect() const
{
    return str_vect_;
}

const void FeatureTabElem::RestoreCountVect(std::ifstream &idx_file, const size_t nb_count)
{
    if (!count_vect_.empty())
    {
        throw std::domain_error("restore non-empty count vector");
    }
    count_vect_.resize(nb_count);
    idx_file.seekg(index_pos_);
    idx_file.read(reinterpret_cast<char *>(&count_vect_[0]), nb_count * sizeof(count_vect_[0]));
}

const void FeatureTabElem::ClearCountVect()
{
    count_vect_.clear();
}
