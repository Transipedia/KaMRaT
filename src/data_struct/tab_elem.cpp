#include <iostream>

#include "tab_elem.hpp"
#include "statistics.hpp"

TabElem::TabElem(std::istringstream &line_conv, std::ofstream &idx_file,
                 float &rep_val,
                 const TabHeader &tab_header)
    : index_pos_(idx_file.is_open() ? static_cast<size_t>(idx_file.tellp()) : 0)
{
    count_vect_.reserve(tab_header.GetNbCount());
    value_vect_.reserve(tab_header.GetNbValue());

    static std::string term;
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    for (unsigned int i(1); line_conv >> term, i < tab_header.GetNbCol(); ++i)
    {
        char nat_ch = tab_header.GetColNatureAt(i);
        if (nat_ch == 'v') // v for values
        {
            value_vect_.emplace_back(std::move(std::stof(term)));
        }
        else if (nat_ch >= 'A' && nat_ch <= 'Z') // A-Z for sample conditions, maximally support 26 conditions
        {
            count_vect_.emplace_back(std::move(std::stof(term)));
        }
        else
        {
            throw std::invalid_argument("unknown column nature code: " + nat_ch);
        }
    }
    static size_t rep_colpos, serial = 0;
    rep_val = (0 == (rep_colpos = tab_header.GetRepColPos()) ? serial++ : value_vect_[tab_header.GetColSerialAt(rep_colpos)]);
    if (line_conv >> term) // should not happen, for debug
    {
        throw std::domain_error("parsing string line to fields failed");
    }
    if (idx_file.is_open())
    {
        idx_file.write(reinterpret_cast<char *>(&count_vect_[0]), count_vect_.size() * sizeof(count_vect_[0])); // first counts
        idx_file.write(reinterpret_cast<char *>(&value_vect_[0]), value_vect_.size() * sizeof(value_vect_[0])); // then values
        std::vector<float>(std::move(count_vect_));                                                             // clear and reallocate the vector
        std::vector<float>(std::move(value_vect_));                                                             // clear and reallocate the vector
    }
    else
    {
        count_vect_.shrink_to_fit();
        value_vect_.shrink_to_fit();
    }
}

const std::vector<float> &TabElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count) const
{
    if (idx_file.is_open())
    {
        count_vect.resize(nb_count);
        idx_file.seekg(index_pos_);
        idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(count_vect[0]));
        return count_vect;
    }
    else
    {
        return count_vect_;
    }
}

const void TabElem::GetVectsAndClear(std::vector<float> &count_vect, std::vector<float> &value_vect,
                                     std::ifstream &idx_file, const size_t nb_count, const size_t nb_value)
{
    if (idx_file.is_open())
    {
        count_vect.clear();
        value_vect.clear();
        count_vect.resize(nb_count);
        value_vect.resize(nb_value);
        idx_file.seekg(index_pos_);
        idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(count_vect[0]));
        idx_file.read(reinterpret_cast<char *>(&value_vect[0]), nb_value * sizeof(value_vect[0]));
    }
    else
    {
        count_vect = std::move(count_vect_);
        value_vect = std::move(value_vect_);
    }
}
