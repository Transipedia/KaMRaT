#include <sstream>
#include <iostream>

#include "count_tab.hpp"
#include "statistics.hpp"

inline void LoadCountFromIndex(std::vector<float> &counts, std::istream &idx_file, const size_t position_on_disk, const size_t nb_samples)
{
    counts.resize(nb_samples);
    idx_file.seekg(position_on_disk);
    idx_file.read(reinterpret_cast<char *>(&counts[0]), nb_samples * sizeof(counts[0]));
}

CountTab::CountTab(const std::string &idx_file_path)
    : idx_file_path_(idx_file_path)
{
}

const bool CountTab::IsCountsInMem() const
{
    return idx_file_path_.empty();
}

const size_t CountTab::GetTableSize() const
{
    return (IsCountsInMem() ? count_tab_.size() : index_pos_.size());
}

const std::string &CountTab::GetIndexPath() const
{
    return idx_file_path_;
}

const bool CountTab::AddRowAsFields(float &row_score, const std::string &line_str, std::ofstream &idx_file)
{
    std::vector<float> value_vect, count_vect;
    std::vector<std::string> str_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, str_vect, line_str);

    str_tab_.emplace_back(str_vect); // non-value columns are alwayes saved in memory
    if (nb_value_ > 0)
    {
        value_tab_.emplace_back(value_vect); // non-count values are always saved in memory
    }

    if (!idx_file.is_open()) // if counts in memory
    {
        count_tab_.emplace_back(count_vect);
    }
    else // if counts on disk
    {
        index_pos_.emplace_back(idx_file.tellp());
        idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(count_vect[0]));
    }

    if (score_pos == 0) // if without score column
    {
        row_score = 0;
        return false;
    }
    else // if with score column
    {
        row_score = value_vect[colserial_vect_[score_pos]];
        return true;
    }
}

const bool CountTab::AddRowAsString(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file)
{
    std::vector<float> value_vect;
    std::vector<std::string> str_vect;
    size_t score_pos = ParseLineStr(value_vect, count_vect, str_vect, line_str);

    index_pos_.emplace_back(idx_file.tellp());
    idx_file << line_str << std::endl;

    if (score_pos == 0)
    {
        row_score = 0;
        return false;
    }
    else
    {
        row_score = value_vect[colserial_vect_[score_pos]];
        return true;
    }
}

const float CountTab::GetValue(const size_t kmer_serial, const size_t valcol_serial) const
{
    return value_tab_.at(kmer_serial).at(valcol_serial);
}

const void CountTab::GetCountVect(std::vector<float> &count_vect, const size_t kmer_serial, std::ifstream &idx_file) const
{
    if (kmer_serial >= str_tab_.size())
    {
        throw std::domain_error("k-mer serial " + std::to_string(kmer_serial) + " larger than count table size " + std::to_string(str_tab_.size()));
    }

    if (idx_file_path_.empty()) // if counts in memory
    {
        count_vect = count_tab_[kmer_serial];
    }
    else if (idx_file.is_open()) // if counts on disk
    {
        LoadCountFromIndex(count_vect, idx_file, index_pos_[kmer_serial], nb_count_);
    }
    else
    {
        throw std::domain_error("searching index on disk but index file not opened");
    }
}

const void CountTab::GetRowString(std::string &row_string, const size_t row_serial, std::ifstream &idx_file, const std::string &as_scorecol) const
{
    size_t disk_pos = index_pos_.at(row_serial);
    std::string str_line;
    if (idx_file.is_open())
    {
        idx_file.seekg(disk_pos);
        std::getline(idx_file, str_line);
    }
    else
    {
        throw std::domain_error("searching index on disk but index file not opened");
    }
    row_string = str_line.substr(0, str_line.find_first_of(" \t") + 1) + as_scorecol + str_line.substr(str_line.find_first_of(" \t"));
}

const float CountTab::AddCountVectIfCoherent(std::vector<float> &sum_count_vect,
                                             const size_t row_serial, const std::vector<float> &ref_count_vect,
                                             const std::string &dist_method, const float dist_thres,
                                             std::ifstream &idx_file) const
{
    float count_dist = 1;
    if (row_serial >= str_tab_.size())
    {
        throw std::domain_error("k-mer serial " + std::to_string(row_serial) + " larger than count table size " + std::to_string(str_tab_.size()));
    }
    if (idx_file_path_.empty()) // if counts in memory
    {
        count_dist = (dist_method == "none" ? 0 : CalcDistance(count_tab_[row_serial], ref_count_vect, dist_method));
        if (dist_method == "none" || count_dist < dist_thres)
        {
            for (size_t col_serial(0); col_serial < nb_count_; ++col_serial)
            {
                sum_count_vect[col_serial] += count_tab_[row_serial][col_serial];
            }
        }
    }
    else if (idx_file.is_open()) // if counts on disk
    {
        std::vector<float> count_vect_x;
        LoadCountFromIndex(count_vect_x, idx_file, index_pos_[row_serial], nb_count_);
        count_dist = (dist_method == "none" ? 0 : CalcDistance(count_vect_x, ref_count_vect, dist_method));
        if (dist_method == "none" || count_dist < dist_thres)
        {
            for (size_t col_serial(0); col_serial < nb_count_; ++col_serial)
            {
                sum_count_vect[col_serial] += count_vect_x[col_serial];
            }
        }
    }
    else
    {
        throw std::domain_error("searching index on disk but index file not opened");
    }
    return count_dist;
}

const float CountTab::CalcCountDistance(const size_t row_serial1, const size_t row_serial2, const std::string &dist_method, std::ifstream &idx_file) const
{
    if (row_serial1 >= str_tab_.size() || row_serial2 >= str_tab_.size())
    {
        std::cerr << row_serial1 << "\t" << row_serial2 << std::endl;
        throw std::domain_error("k-mer serial larger than count table size " + std::to_string(str_tab_.size()));
    }

    if (idx_file_path_.empty()) // if counts in memory
    {
        return CalcDistance(count_tab_[row_serial1], count_tab_[row_serial2], dist_method);
    }
    else if (idx_file.is_open()) // if counts on disk
    {
        std::vector<float> count_vect1, count_vect2;
        LoadCountFromIndex(count_vect1, idx_file, index_pos_[row_serial1], nb_count_);
        LoadCountFromIndex(count_vect2, idx_file, index_pos_[row_serial2], nb_count_);
        return CalcDistance(count_vect1, count_vect2, dist_method);
    }
    else
    {
        throw std::domain_error("searching index on disk but index file not opened");
    }
}
