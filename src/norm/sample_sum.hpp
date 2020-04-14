#ifndef KAMRAT_NORM_SAMPLESUM_HPP
#define KAMRAT_NORM_SAMPLESUM_HPP

#include <string>

class SampleSum
{
public:
    SampleSum(const std::string &sample_name);
    void AddCount(double count);
    double GetCount() const;
    void Print(std::ostream &out_s) const;

private:
    std::string sample_name_;
    double sum_count_;
};

#endif //KAMRAT_NORM_SAMPLESUM_HPP