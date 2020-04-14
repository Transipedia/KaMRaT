#include <iostream>

#include "sample_sum.hpp"

SampleSum::SampleSum(const std::string &sample_name)
    : sample_name_(sample_name), sum_count_(0)
{
}

void SampleSum::AddCount(const double count)
{
    sum_count_ += count;
}

void SampleSum::Print(std::ostream &out_s) const
{
    out_s << sample_name_ << "\t" << sum_count_ << std::endl;
}

double SampleSum::GetCount() const
{
    return sum_count_;
}
