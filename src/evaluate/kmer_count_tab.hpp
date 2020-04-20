#ifndef KAMRAT_EVALUATE_KMERCOUNTTAB_HPP
#define KAMRAT_EVALUATE_KMERCOUNTTAB_HPP

#include <map>
#include <vector>

class KMerCountTab
{
public:
    const bool AddKMerCountInMem(uint64_t kmer_code, const std::vector<float> &kmer_counts); // return whether adding is succeeded
    const bool CheckKMerCodeExist(uint64_t kmer_code) const;
    const std::vector<float> &GetCountByKmer(uint64_t kmer_code) const;
private:
    std::map<uint64_t, size_t> kmer_count_dict_;
    std::vector<std::vector<float>> count_tab_;
};

#endif //KAMRAT_EVALUATE_KMERCOUNTTAB_HPP