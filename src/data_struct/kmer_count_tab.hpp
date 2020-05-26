#ifndef KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP
#define KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP

#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <map>

class KMerCountTab
{
public:
    KMerCountTab(const std::string &mode);
    const std::string GetMode() const;
    //----- header info -----//
    const void MakeColumnInfo(const std::string &header_line, const std::string &sample_info_path, const std::string &score_colname);
    const std::string GetColName(size_t i_col) const;
    const char GetColNature(size_t i_col) const;
    const size_t GetColSerial(size_t i_col) const;
    const int GetSmpLabel(size_t i_col) const;
    const size_t GetNbCondition() const;
    const size_t GetNbValue() const;
    const size_t GetNbCount() const;
    const size_t GetNbColumn() const;
    //----- value-count info -----//
    const float AddKMerCountInMem(const std::string &line_str);
    const float AddKMerIndexOnDsk(const std::string &line_str, std::ofstream &index_file);
    const float GetValue(size_t kmer_serial, size_t valcol_serial) const;
    const void GetCountInMem(std::vector<float> &count_vect, size_t kmer_serial) const;
    const void GetCountOnDsk(std::vector<float> &count_vect, size_t kmer_serial, std::ifstream &index_file) const;
    const size_t GetAvgCountInMem(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set) const;
    const size_t GetAvgCountOnDsk(std::vector<float> &count_avg_vect, const std::set<size_t> kmer_serial_set, std::ifstream &index_file) const;

private:
    const std::string mode_; // inMem or onDsk
    //----- header info -----//
    size_t nb_cond_, nb_value_, nb_count_;  // number of conditions and value/count columns
    std::vector<std::string> colname_vect_; // column name vector parsed from the header line
    std::vector<char> colnature_vect_;      // column nature vector, f=feature s=sample v=value +=score_value
    std::vector<size_t> colserial_vect_;    // column serial vector, serial number is assigned SEPARATELY to count or to value
    std::vector<size_t> smplabel_vect_;     // sample label vector parsed from the header line and sample-info file
    //----- value-count info -----//
    std::vector<std::vector<float>> value_tab_; // k-mer value vector ------ used both in inMem or onDsk
    std::vector<std::vector<float>> count_tab_; // k-mer count vector ------ only used in inMem
    std::vector<size_t> index_pos_;             // indexed position for k-mer count ------ only used in onDsk
};

using code2serial_t = std::map<uint64_t, size_t>; // external dictionary to link k-mer with row serial number

#endif //KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP