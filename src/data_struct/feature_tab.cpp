#include <sstream>
#include <iostream>

#include "feature_tab.hpp"
#include "statistics.hpp"

FeatureTab::FeatureTab(const std::string &idx_file_path,
                       const std::string &sample_info_path,
                       const std::unordered_set<std::string> &preserved_cond_tags = std::unordered_set<std::string>())
    : idx_file_path_(idx_file_path),
      nb_row_(0),
      feature_tab_header_(sample_info_path, preserved_cond_tags)
{
    std::cout << "feature tab initialized: "
              << feature_tab_header_.GetNbCol()
              << feature_tab_header_.GetNbCondition()
              << feature_tab_header_.GetNbValue()
              << feature_tab_header_.GetNbCount()
              << feature_tab_header_.GetNbStr() << std::endl;
}

const std::string &FeatureTab::GetIndexPath() const
{
    return idx_file_path_;
}

const size_t FeatureTab::GetTableSize() const
{
    return nb_row_;
}

const void FeatureTab::MakeTable(const std::string &table_path, const std::string &repval_name)
{
    
}

// const void FeatureTab::MakeHeaderRow(const std::string &header_line_str, const std::string &rep_colname)
// {
//     repval_colpos_ = feature_tab_header_.MakeColumnInfo(header_line_str, rep_colname);
// }

// const bool CountTab::AddRowAsFields(float &row_score, const std::string &line_str, std::ofstream &idx_file)
// {
//     feature_tab_.emplace_back(tab_header_.GetNbValue(), tab_header_.GetNbCount(), tab_header_.GetNbStr(), std::istringstream(line_str), tab_header_.GetColNatureVect());

//     if (nb_value_ > 0)
//     {
//         value_tab_.emplace_back(std::move(value_vect)); // non-count values are always saved in memory
//     }
//     if (nb_str_ > 0)
//     {
//         str_tab_.emplace_back(std::move(str_vect)); // non-value columns are alwayes saved in memory
//     }

//     if (!idx_file.is_open()) // if counts in memory
//     {
//         count_tab_.emplace_back(std::move(count_vect));
//     }
//     else // if counts on disk
//     {
//         index_pos_.emplace_back(idx_file.tellp());
//         idx_file.write(reinterpret_cast<char *>(&count_vect[0]), count_vect.size() * sizeof(count_vect[0]));
//     }

//     value_vect.clear();
//     count_vect.clear();
//     str_vect.clear();

//     if (score_pos == 0) // if without score column
//     {
//         row_score = 0;
//         return false;
//     }
//     else // if with score column
//     {
//         row_score = value_tab_.back()[colserial_vect_[score_pos]];
//         return true;
//     }
// }

// const bool CountTab::AddRowAsString(float &row_score, std::vector<float> &count_vect, const std::string &line_str, std::ofstream &idx_file)
// {
//     std::vector<float> value_vect;
//     std::vector<std::string> str_vect;
//     size_t score_pos = ParseLineStr(value_vect, count_vect, str_vect, line_str);

//     index_pos_.emplace_back(idx_file.tellp());
//     idx_file << std::move(line_str) << std::endl;

//     if (score_pos == 0)
//     {
//         row_score = 0;
//         return false;
//     }
//     else
//     {
//         row_score = value_vect[colserial_vect_[score_pos]];
//         return true;
//     }
// }

// const float CountTab::GetValue(const size_t kmer_serial, const size_t valcol_serial) const
// {
//     return value_tab_.at(kmer_serial).at(valcol_serial);
// }

// const void CountTab::GetCountVect(std::vector<float> &count_vect, const size_t kmer_serial, std::ifstream &idx_file) const
// {
//     if (kmer_serial >= (IsCountsInMem() ? count_tab_.size() : index_pos_.size()))
//     {
//         throw std::domain_error("k-mer serial " + std::to_string(kmer_serial) + " larger than count table size " + std::to_string(str_tab_.size()));
//     }

//     if (idx_file_path_.empty()) // if counts in memory
//     {
//         count_vect = count_tab_[kmer_serial];
//     }
//     else if (idx_file.is_open()) // if counts on disk
//     {
//         LoadCountFromIndex(count_vect, idx_file, index_pos_[kmer_serial], nb_count_);
//     }
//     else
//     {
//         throw std::domain_error("searching index on disk but index file not opened");
//     }
// }

// const void CountTab::GetRowString(std::string &row_string, const size_t row_serial, std::ifstream &idx_file, const std::string &as_scorecol) const
// {
//     size_t disk_pos = index_pos_.at(row_serial);
//     std::string str_line;
//     if (idx_file.is_open())
//     {
//         idx_file.seekg(disk_pos);
//         std::getline(idx_file, str_line);
//     }
//     else
//     {
//         throw std::domain_error("searching index on disk but index file not opened");
//     }
//     row_string = str_line.substr(0, str_line.find_first_of(" \t") + 1) + as_scorecol + str_line.substr(str_line.find_first_of(" \t"));
// }

// const void CountTab::EstimateMeanCountVect(std::vector<float> &avg_count_vect, const std::vector<size_t> &row_serial_vect, std::ifstream &idx_file) const
// {
//     avg_count_vect.assign(nb_count_, 0);
//     for (size_t rs : row_serial_vect)
//     {
//         if (rs >= (IsCountsInMem() ? count_tab_.size() : index_pos_.size())) // should not happen
//         {
//             throw std::domain_error("k-mer serial " + std::to_string(rs) + " larger than count table size " + std::to_string(str_tab_.size()));
//         }
//         if (idx_file_path_.empty()) // if counts in memory
//         {
//             for (size_t col_serial(0); col_serial < nb_count_; ++col_serial)
//             {
//                 avg_count_vect[col_serial] += count_tab_[rs][col_serial];
//             }
//         }
//         else if (idx_file.is_open()) // if counts on disk
//         {
//             static std::vector<float> count_vect_x;
//             LoadCountFromIndex(count_vect_x, idx_file, index_pos_[rs], nb_count_);
//             for (size_t col_serial(0); col_serial < nb_count_; ++col_serial)
//             {
//                 avg_count_vect[col_serial] += count_vect_x[col_serial];
//             }
//         }
//         else
//         {
//             throw std::domain_error("searching index on disk but index file not opened");
//         }
//     }
//     for (size_t i(0); i < avg_count_vect.size(); ++i)
//     {
//         avg_count_vect[i] /= row_serial_vect.size();
//     }
// }

// const float CountTab::CalcCountDistance(const size_t row_serial1, const size_t row_serial2, const std::string &dist_method, std::ifstream &idx_file) const
// {
//     size_t tab_size = (IsCountsInMem() ? count_tab_.size() : index_pos_.size());
//     if (row_serial1 >= tab_size || row_serial2 >= tab_size)
//     {
//         std::cerr << row_serial1 << "\t" << row_serial2 << std::endl;
//         throw std::domain_error("k-mer serial larger than count table size " + std::to_string(str_tab_.size()));
//     }

//     if (idx_file_path_.empty()) // if counts in memory
//     {
//         return CalcDistance(count_tab_[row_serial1], count_tab_[row_serial2], dist_method);
//     }
//     else if (idx_file.is_open()) // if counts on disk
//     {
//         static std::vector<float> count_vect1, count_vect2;
//         LoadCountFromIndex(count_vect1, idx_file, index_pos_[row_serial1], nb_count_);
//         LoadCountFromIndex(count_vect2, idx_file, index_pos_[row_serial2], nb_count_);
//         return CalcDistance(count_vect1, count_vect2, dist_method);
//     }
//     else
//     {
//         throw std::domain_error("searching index on disk but index file not opened");
//     }
// }

// const void CountTab::ShrinkTab()
// {
//     value_tab_.shrink_to_fit();
//     count_tab_.shrink_to_fit();
//     str_tab_.shrink_to_fit();
//     index_pos_.shrink_to_fit();
// }
