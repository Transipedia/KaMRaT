#ifndef KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP
#define KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP

#include <map>
#include <vector>
#include <string>

template <typename countT>
class KMerCountTab
{
public:
    const bool AddKMerCountInMem(uint64_t kmer_code, const std::vector<countT> &kmer_counts); // return whether adding is succeeded
    const bool AddKMerIndexOnDsk(uint64_t kmer_code, size_t dsk_pos); // return whether adding is succeeded
    const bool IsKMerExist(uint64_t kmer_code) const;
    const bool GetCountInMem(std::vector<countT> &count, uint64_t kmer_code) const;
    const bool GetIndexOnDsk(size_t &index_val, uint64_t kmer_code) const;
private:
    std::map<uint64_t, size_t> kmer_count_dict_;
    std::vector<std::vector<countT>> count_tab_;
};

#endif //KAMRAT_DATASTRUCT_KMERCOUNTTAB_HPP