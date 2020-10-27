#include <iostream>

#include "feature_tab_elem.hpp"
#include "statistics.hpp"

FeatureTabElem::FeatureTabElem(std::istringstream &line_conv, std::ofstream &idx_file, const FeatureTabHeader &tab_header)
    : index_pos_(idx_file.is_open() ? static_cast<size_t>(idx_file.tellp()) : 0)
{
    value_vect_.reserve(tab_header.GetNbValue());
    count_vect_.reserve(tab_header.GetNbCount());
    str_vect_.reserve(tab_header.GetNbStr());

    static std::string term;
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    char nat_ch;
    for (unsigned int i(1); line_conv >> term, i < tab_header.GetNbCol(); ++i)
    {
        nat_ch = tab_header.GetColNatureAt(i);
        if (nat_ch == 'v') // v for values
        {
            value_vect_.emplace_back(std::move(std::stof(term)));
        }
        else if (nat_ch >= 'A' && nat_ch <= 'Z') // A-Z for sample conditions, maximally support 26 conditions
        {
            count_vect_.emplace_back(std::move(std::stof(term)));
        }
        else if (nat_ch == 's') // s for strings
        {
            str_vect_.emplace_back(std::move(term));
        }
        else
        {
            throw std::invalid_argument("unknown column nature code: " + nat_ch);
        }
    }
    if (line_conv >> term) // should not happen, for debug
    {
        throw std::domain_error("parsing string line to fields failed");
    }
    // for (unsigned int i(1); line_conv >> term, i < tab_header.GetNbCol(); ++i)
    // {
    //     std::cout << "\t" << count_vect_[tab_header.GetColSerialAt(i)];
    // }
    // std::cout << std::endl;
    if (idx_file.is_open())
    {
        idx_file.write(reinterpret_cast<char *>(&value_vect_[0]), value_vect_.size() * sizeof(value_vect_[0])); // first values
        idx_file.write(reinterpret_cast<char *>(&count_vect_[0]), count_vect_.size() * sizeof(count_vect_[0])); // then counts
        for (size_t i(0); i < tab_header.GetNbStr(); ++i)                                                       // finally strings
        {
            idx_file << str_vect_[i];
        }
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

const void FeatureTabElem::RestoreRow(std::ifstream &idx_file, const size_t nb_value, const size_t nb_count, const size_t nb_str)
{
    if (!value_vect_.empty() || !count_vect_.empty() || !str_vect_.empty()) // should not happen, for debug
    {
        throw std::domain_error("restore non-empty count/value/string vector");
    }
    value_vect_.resize(nb_value);
    count_vect_.resize(nb_count);
    str_vect_.resize(nb_str);
    idx_file.seekg(index_pos_);
    idx_file.read(reinterpret_cast<char *>(&value_vect_[0]), nb_value * sizeof(value_vect_[0])); // first values
    idx_file.read(reinterpret_cast<char *>(&count_vect_[0]), nb_count * sizeof(count_vect_[0])); // then counts
    for (size_t i(0); i < nb_str; ++i)                                                           // finally strings
    {
        idx_file >> str_vect_[i];
    }
}

const void FeatureTabElem::ClearRow()
{
    std::vector<float>(std::move(value_vect_));
    std::vector<float>(std::move(count_vect_));
    std::vector<std::string>(std::move(str_vect_));
}
