#ifndef KAMRAT_UTILS_SAMPLEINFO_HPP
#define KAMRAT_UTILS_SAMPLEINFO_HPP

#include <iostream>
#include <unordered_map>
#include <string>

class SampleInfo
{
public:
    SampleInfo();
    void AddSampleInfo(const std::string &sample_info_line);
    const bool IsEmpty() const;
    const size_t GetNbCondition() const;
    const size_t GetNbSample() const;
    const int GetLabel(const std::string &sample_name) const;
    const bool IsSample(const std::string &col_name) const;
    void PrintSampleNames(std::ostream &out_s, const std::string &leading_str) const;

private:
    size_t nb_condition_;
    std::unordered_map<std::string, int> sample_label_;
    std::unordered_map<std::string, int> condition_label_;
};

#endif //KAMRAT_UTILS_SAMPLEINFO_HPP