#include <iostream>

#include "sample_info.hpp"

SampleInfo::SampleInfo(const std::string &sample_name)
    : sample_name_(sample_name), sample_count_(0), sample_label_(0) // TODO: make it coherent with kamratReduce
{
}

const void SampleInfo::AddCount(const double count)
{
    sample_count_ += count;
}

const double SampleInfo::GetCount() const
{
    return sample_count_;
}