#include <fstream>
#include <sstream>

#include "tab_header.hpp"

TabHeader::TabHeader(const std::string &sample_info_path, std::unordered_set<std::string> &&preserved_cond_tags)
    : nb_value_(0), nb_count_(0), nb_str_(0), nb_col_(0), nb_condi_(0), rep_colpos_(0)
{
    if (!sample_info_path.empty())
    {
        std::ifstream sample_info_file(sample_info_path);
        if (!sample_info_file.is_open())
        {
            throw std::domain_error("could not open sample-info file: " + sample_info_path);
        }
        std::string sample_info_line, sample, condition;
        std::istringstream conv;
        while (std::getline(sample_info_file, sample_info_line))
        {
            conv.str(sample_info_line);
            conv >> sample >> condition;
            if (conv.fail())
            {
                condition = ""; // empty string indicating the unique condition
            }
            if (preserved_cond_tags.find(condition) != preserved_cond_tags.cend())
            {
                throw std::domain_error("condition name (" + condition + ") is preserved by KaMRaT, please change another name for this condition");
            }
            auto ins = condi2lab_.insert({condition, 'A' + nb_condi_}); // try to add the condition into the condition dictionary
            if (ins.second)                                             // if a new condition is added to the dictionary
            {
                ++nb_condi_;
            }
            if (!smp2lab_.insert({sample, ins.first->second}).second) // associate the sample with its condition label
            {
                throw std::domain_error("sample-info file has duplicated sample name: " + sample);
            }
            conv.clear();
        }
        sample_info_file.close();
        if (nb_condi_ != condi2lab_.size()) // should not happen, for debug
        {
            throw std::domain_error("number of condition not consistent");
        }
    }
}

const void TabHeader::MakeColumnInfo(const std::string &header_line, const std::string &rep_colname)
{
    if (header_line.empty()) // quit if header line is empty
    {
        throw std::domain_error("cannot parse column information with an empty header line");
    }
    std::istringstream conv(header_line);
    std::string term;
    // the first column is supposed to be feature string //
    conv >> term;
    colname_vect_.emplace_back(term);
    colserial_vect_.emplace_back(nb_col_++); // place hoder: first column always as string
    colnature_vect_.emplace_back('s'); // place hoder: first column always as string
    // following columns //
    if (condi2lab_.empty()) // if sample info file NOT provided, then all next columns are samples
    {
        while (conv >> term)
        {
            colname_vect_.emplace_back(term);
            colnature_vect_.emplace_back('A');
            colserial_vect_.emplace_back(nb_count_++);
            ++nb_col_;
        }
    }
    else // sample info file provided => next columns may contain non-sample values
    {
        while (conv >> term)
        {
            colname_vect_.emplace_back(term);
            auto iter = smp2lab_.find(term);
            if (iter != smp2lab_.cend())
            {
                colnature_vect_.emplace_back(iter->second);
                colserial_vect_.emplace_back(nb_count_++);
            }
            else
            {
                if (term == rep_colname)
                {
                    rep_colpos_ = nb_col_;
                }
                colnature_vect_.emplace_back('v');
                colserial_vect_.emplace_back(nb_value_++);
            }
            ++nb_col_;
        }
    }
    if (nb_count_ == 0)
    {
        throw std::domain_error("no sample column found");
    }
    if (nb_col_ != colname_vect_.size() || nb_col_ != colserial_vect_.size() || nb_col_ != colnature_vect_.size()) // should not happen, for debug
    {
        throw std::domain_error("colname, colserial, colnature sizes not coherent");
    }
}

const size_t TabHeader::GetNbValue() const
{
    return nb_value_;
}

const size_t TabHeader::GetNbCount() const
{
    return nb_count_;
}

const size_t TabHeader::GetNbStr() const
{
    return nb_str_;
}

const size_t TabHeader::GetNbCol() const
{
    return nb_col_;
}

const size_t TabHeader::GetNbCondition() const
{
    return nb_condi_;
}

const size_t TabHeader::GetRepColPos() const
{
    return rep_colpos_;
}

const char TabHeader::GetConditionLabel(const std::string &condi) const
{
    auto iter = condi2lab_.find(condi);
    return (iter == condi2lab_.cend() ? iter->second : '\0');
}

const float TabHeader::ParseRowStr(std::vector<float> &count_vect, std::vector<float> &value_vect, std::istringstream &line_conv) const
{
    count_vect.reserve(nb_count_);
    value_vect.reserve(nb_value_);

    static std::string term;
    line_conv >> term; // skip the first column which is k-mer/tag/contig/sequence
    for (unsigned int i(1); i < nb_col_ && (line_conv >> term); ++i)
    {
        char nat_ch = colnature_vect_[i];
        if (nat_ch == 'v') // v for values
        {
            value_vect.emplace_back(std::stof(std::move(term)));
        }
        else if (nat_ch >= 'A' && nat_ch <= 'Z') // A-Z for sample conditions, maximally support 26 conditions
        {
            count_vect.emplace_back(std::stof(std::move(term)));
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
    return (0 == rep_colpos_ ? 0 : value_vect[colserial_vect_[rep_colpos_]]);
}

const std::string &TabHeader::GetColNameAt(const size_t i) const
{
    if (i > nb_col_)
    {
        throw std::domain_error("getting column name with an index exceeds column number");
    }
    return colname_vect_[i];
}

const char TabHeader::GetColNatureAt(const size_t i) const
{
    if (i > nb_col_)
    {
        throw std::domain_error("getting column nature with an index exceeds column number");
    }
    return colnature_vect_[i];
}

const size_t TabHeader::GetColSerialAt(const size_t i) const
{
    if (i > nb_col_)
    {
        throw std::domain_error("getting column serial with an index exceeds column number");
    }
    return colserial_vect_[i];
}

const bool TabHeader::IsSample(size_t i_col) const
{
    if (i_col > nb_col_)
    {
        throw std::domain_error("column number in search exceeds table's column number");
    }
    return (colnature_vect_[i_col] >= 'A' && colnature_vect_[i_col] <= 'Z');
}

const void TabHeader::ParseSmpLabels(std::vector<size_t> &smp_labels) const
{
    if (!smp_labels.empty())
    {
        throw std::domain_error("sample label vector not empty before adding sample labels");
    }
    for (size_t i(0); i < colnature_vect_.size(); ++i)
    {
        if (colnature_vect_[i] >= 'A' && colnature_vect_[i] <= 'Z')
        {
            smp_labels.emplace_back(colnature_vect_[i] - 'A');
        }
    }
}
