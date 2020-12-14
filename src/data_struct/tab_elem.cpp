#include <iostream>

#include "tab_elem.hpp"

TabElem::TabElem(std::istringstream &line_conv, std::ofstream &idx_file, const TabHeader &tab_header)
{
    idx_pos_ = (idx_file.is_open() ? static_cast<size_t>(idx_file.tellp()) : 0);
    value_ = tab_header.ParseRowStr(count_vect_, non_count_str_, line_conv);
    if (idx_file.is_open())
    {
        idx_file.write(reinterpret_cast<char *>(&count_vect_[0]), count_vect_.size() * sizeof(count_vect_[0])); // first counts columns
        idx_file << non_count_str_ << std::endl;                                                                // then non-count columns
        std::vector<float>().swap(count_vect_);                                                                 // clear and reallocate the vector
        std::string().swap(non_count_str_);                                                                     // clear and reallocate the vector
        // std::cout << count_vect_.capacity() << "\t" << non_count_str_.capacity() << std::endl;
    }
    else
    {
        count_vect_.shrink_to_fit();
        non_count_str_.shrink_to_fit();
    }
}

TabElem::TabElem(std::istringstream &line_conv, std::vector<float> &count_vect, std::ofstream &idx_file, const TabHeader &tab_header)
{
    idx_pos_ = (idx_file.is_open() ? static_cast<size_t>(idx_file.tellp()) : 0);
    value_ = tab_header.ParseRowStr(count_vect, non_count_str_, line_conv);
    if (idx_file.is_open())
    {
        idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(count_vect[0])); // first counts columns
        idx_file << non_count_str_ << std::endl;                                                             // then non-count columns
        std::string().swap(non_count_str_);                                                                  // clear and reallocate the vector
    }
    else // should not happen, for debug
    {
        throw std::domain_error("Forcing indexing without index file opened");
    }
}

const std::string &TabElem::MakeOutputRowStr(std::string &row_str, std::ifstream &idx_file, const size_t nb_count)
{
    row_str.clear();
    if (idx_file.is_open())
    {
        count_vect_.resize(nb_count);
        idx_file.seekg(idx_pos_);
        idx_file.read(reinterpret_cast<char *>(&count_vect_[0]), nb_count * sizeof(count_vect_[0]));
        std::getline(idx_file, row_str);
    }
    for (size_t i(0); i < nb_count; ++i)
    {
        row_str += ("\t" + std::to_string(count_vect_[i]));
    }
    return row_str;
}

const std::string &TabElem::MakeOutputRowStr(std::string &row_str, const std::vector<float> &count_vect, std::ifstream &idx_file)
{
    row_str.clear();
    const size_t nb_count = count_vect.size();
    GetNonCountStr(row_str, idx_file, nb_count);
    for (size_t i(0); i < nb_count; ++i)
    {
        row_str += ("\t" + std::to_string(count_vect[i]));
    }
    return row_str;
}

const float TabElem::GetValue() const
{
    return value_;
}

const std::vector<float> &TabElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count) const
{
    count_vect.clear();
    if (idx_file.is_open())
    {
        count_vect.resize(nb_count);
        idx_file.seekg(idx_pos_);
        idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(count_vect[0]));
    }
    else
    {
        count_vect = count_vect_;
    }
    return count_vect;
}

const std::string &TabElem::GetNonCountStr(std::string &non_count_str, std::ifstream &idx_file, const size_t nb_count) const
{
    non_count_str.clear();
    if (idx_file.is_open())
    {
        idx_file.seekg(idx_pos_ + nb_count * sizeof(decltype(count_vect_)::value_type));
        std::getline(idx_file, non_count_str);
    }
    else
    {
        non_count_str = non_count_str_;
    }
    return non_count_str;
}

const size_t TabElem::GetIdxPos() const
{
    return idx_pos_;
}