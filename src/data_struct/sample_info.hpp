#ifndef KAMRAT_DATASTRUCT_SAMPLEINFO_HPP
#define KAMRAT_DATASTRUCT_SAMPLEINFO_HPP

#include <string>
#include <vector>

class SampleInfo
{
public:
    SampleInfo(const std::string &sample_name);
    const void AddCount(double count);
    const double GetCount() const;

private:
    std::string sample_name_;
    double sample_count_;
    size_t sample_label_;
};

using sampleInfoVect_t = std::vector<SampleInfo>;

#endif //KAMRAT_DATASTRUCT_SAMPLEINFO_HPP