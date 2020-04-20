#include "kmer_count_tab.hpp"

const bool KMerCountTab::AddKMerCountInMem(const uint64_t kmer_code, const std::vector<float> &kmer_counts)
{
    bool is_succeeded = kmer_count_dict_.insert({kmer_code, count_tab_.size()}).second;
    count_tab_.emplace_back(kmer_counts);
    return is_succeeded;
}

const bool KMerCountTab::CheckKMerCodeExist(const uint64_t kmer_code) const
{
    return (kmer_count_dict_.find(kmer_code) != kmer_count_dict_.cend());
}

// Should assure before calling that the k-mer code does exist in k-mer count dictionary ! //
const std::vector<float> &KMerCountTab::GetCountByKmer(const uint64_t kmer_code) const
{
    auto iter = kmer_count_dict_.find(kmer_code);
    return count_tab_.at(iter->second);
}