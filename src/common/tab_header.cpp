#include "tab_header.hpp"

TabHeader::TabHeader()
    : rep_colpos_(0)
{
}

TabHeader::TabHeader(const std::string &smp_info_path)
    : TabHeader()
{
    std::ifstream smp_info_file(smp_info_path);
    if (!smp_info_file.is_open())
    {
        throw std::domain_error("error open sample-info file: " + smp_info_path);
    }
    std::string line, smp_name, condi_name;
    std::istringstream conv;

    size_t condi_serial;
    while (std::getline(smp_info_file, line))
    {
        conv.str(line);
        conv >> smp_name >> condi_name;
        if (conv.fail())
        {
            condi_name = ""; // empty string indicating the unique condition
        }
        for (condi_serial = 0; condi_serial < condi_name_vect_.size(); ++condi_serial)
        {
            if (condi_name_vect_[condi_serial] == condi_name)
            {
                break;
            }
        }
        if (condi_serial == condi_name_vect_.size())
        {
            condi_name_vect_.emplace_back(condi_name);
        }
        if (!smp2lab_.insert({smp_name, condi_serial + 1}).second) // sample labels start from 1, but condition name vector index starts from 0
        {
            throw std::domain_error("sample-info file has duplicated sample name: " + smp_name);
        }
        conv.clear();
    }
    smp_info_file.close();
}

const void TabHeader::MakeColumns(std::istringstream &line_conv, const std::string &rep_colname)
{
    std::string term;
    line_conv >> term;                                       // first column is supposed to be feature string
    col_namenat_vect_.emplace_back(std::make_pair(term, 0)); // not a sample

    if (smp2lab_.empty()) // following columns without given sample-info file
    {
        while (line_conv >> term)
        {
            col_namenat_vect_.emplace_back(std::make_pair(term, 1)); // sample without given condition label
        }
    }
    else // following columns with given sample-info file
    {
        while (line_conv >> term)
        {
            if (rep_colname == term)
            {
                rep_colpos_ = col_namenat_vect_.size(); // calculate before emplace_back(), no need minus 1
                col_namenat_vect_.emplace_back(std::make_pair(term, 0));
            }
            else
            {
                std::unordered_map<std::string, size_t>::const_iterator &&iter = smp2lab_.find(term);
                size_t col_label = (iter == smp2lab_.cend() ? 0 : iter->second);
                col_namenat_vect_.emplace_back(std::make_pair(term, col_label));
            }
        }
    }
    if (!rep_colname.empty() && rep_colpos_ == 0)
    {
        throw std::domain_error("could not find column name: " + rep_colname);
    }
}

const size_t TabHeader::GetNbCondition() const
{
    return (smp2lab_.empty() ? 1 : condi_name_vect_.size()); // if sample-info file not provided, set condition number to 1
}

const size_t TabHeader::GetNbCol() const
{
    return col_namenat_vect_.size();
}

const size_t TabHeader::GetNbCount() const
{
    return (smp2lab_.empty() ? (GetNbCol() - 1) : smp2lab_.size()); // if sample-info file not provided, all except the first column are samples
}

const size_t TabHeader::GetRepColPos() const
{
    return rep_colpos_;
}

const std::string &TabHeader::GetColNameAt(const size_t i_col) const
{
    return col_namenat_vect_[i_col].first;
}

const size_t TabHeader::GetColNatAt(const size_t i_col) const
{
    return col_namenat_vect_[i_col].second;
}

const std::string &TabHeader::GetColCondiAt(const size_t i_col) const
{
    if (!IsColCount(i_col))
    {
        throw std::domain_error("not a sample column");
    }
    else
    {
        return condi_name_vect_[col_namenat_vect_[i_col].second - 1]; // sample labels start from 1, but condition name vector index starts from 0
    }
}

const bool TabHeader::IsColCount(const size_t i_col) const
{
    return (GetColNatAt(i_col) > 0);
}

const std::vector<size_t> &TabHeader::GetSampleLabelVect(std::vector<size_t> &label_vect) const
{
    label_vect.clear();
    for (size_t i_col(1); i_col < GetNbCol(); ++i_col)
    {
        size_t colnat = GetColNatAt(i_col);
        if (colnat > 0)
        {
            label_vect.push_back(colnat - 1);
        }
    }
    return label_vect;
}

const std::vector<std::string> &TabHeader::GetCondiNameVect() const
{
    return condi_name_vect_;
}

const double TabHeader::ParseRowStr(std::vector<float> &count_vect, std::string &non_count_str, std::istringstream &line_conv) const
{
    static std::string term;
    term.clear();
    size_t nb_count = GetNbCount();
    count_vect.clear();
    non_count_str.clear();

    double rep_val(0);
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    non_count_str += term;
    for (size_t i(1); (i < GetNbCol()) && (line_conv >> term); ++i)
    {
        if (IsColCount(i))
        {
            count_vect.push_back(std::stof(term));
        }
        else
        {
            non_count_str += ("\t" + term);
            if (i == rep_colpos_)
            {
                rep_val = std::stod(term);
            }
        }
    }
    if (count_vect.size() != nb_count)
    {
        throw std::domain_error("parsing string line to fields failed: " +
                                std::to_string(count_vect.size()) + " vs " + std::to_string(nb_count));
    }
    return rep_val;
}
