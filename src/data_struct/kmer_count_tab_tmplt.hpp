#ifndef KAMRAT_DATASTRUCT_KMERCOUNTTABTMPLT_HPP
#define KAMRAT_DATASTRUCT_KMERCOUNTTABTMPLT_HPP

#include <map>
#include <vector>
#include <string>

template <typename countT>
class KMerCountTab
{
public:
    const bool AddKMerCountInMem(uint64_t kmer_code, const std::vector<countT> &kmer_counts); // return whether adding is succeeded
    const bool AddKMerIndexOnDsk(uint64_t kmer_code, size_t dsk_pos);                         // return whether adding is succeeded
    const bool IsKMerExist(uint64_t kmer_code) const;
    const bool GetCountInMem(std::vector<countT> &count, uint64_t kmer_code) const;
    const bool GetIndexOnDsk(size_t &index_val, uint64_t kmer_code) const;

private:
    std::map<uint64_t, size_t> kmer_count_dict_;
    std::vector<std::vector<countT>> count_tab_;
};


template <typename countT>
const bool KMerCountTab<countT>::AddKMerCountInMem(const uint64_t kmer_code, const std::vector<countT> &kmer_counts)
/* -------------------------------------------------------------------------------------------- *\
    Args:   1. k-mer unique code
            2. k-mer count vector for insertion
    Value:  bool indicating whether the insertion was succeeded
    Func:   1. check whether the count vector to insert is coherent with the table
            2. insert the k-mer count vector to k-mer count table
            3. associate the k-mer unique code with its inserted row number in the count table
\* -------------------------------------------------------------------------------------------- */
{
    if (count_tab_.empty() || count_tab_.back().size() == kmer_counts.size())
    {
        bool is_succeeded = kmer_count_dict_.insert({kmer_code, count_tab_.size()}).second;
        if (is_succeeded)
        {
            count_tab_.emplace_back(kmer_counts);
        }
        return is_succeeded;
    }
    else
    {
        return false;
    }
}

template <typename countT>
const bool KMerCountTab<countT>::AddKMerIndexOnDsk(const uint64_t kmer_code, const size_t dsk_pos)
/* -------------------------------------------------------------------------------------- *\
    Args:   1. k-mer unique code
            2. indexed position of k-mer count on disk
    Value:  bool indicating whether the insertion was succeeded
    Func:   associate the k-mer unique code with indexed position of k-mer count on disk
\* -------------------------------------------------------------------------------------- */
{
    bool is_succeeded = kmer_count_dict_.insert({kmer_code, dsk_pos}).second;
    return is_succeeded;
}

template <typename countT>
const bool KMerCountTab<countT>::IsKMerExist(const uint64_t kmer_code) const
/* ---------------------------------------------------------------------------- *\
    Arg:    k-mer unique code
    Value   bool indicating whether the k-mer exists in k-mer count dictionary
    Func:   check if the k-mer has been registered in the table
\* ---------------------------------------------------------------------------- */
{
    return (kmer_count_dict_.find(kmer_code) != kmer_count_dict_.cend());
}

template <typename countT>
const bool KMerCountTab<countT>::GetCountInMem(std::vector<countT> &count, const uint64_t kmer_code) const
/* ------------------------------------------------------------------------------------ *\
    Arg:    k-mer unique code
    Value:  a pair object indicating
                1. k-mer counts vectors for storing target counts in query (as refArg)
                2. bool for whether the query was succeeded
    Func:   return the k-mer counts stored in memory by k-mer unique code
\* ------------------------------------------------------------------------------------ */
{
    auto iter = kmer_count_dict_.find(kmer_code);
    if (iter != kmer_count_dict_.cend())
    {
        count = count_tab_.at(iter->second);
        return true;
    }
    else
    {
        return false;
    }
}

template <typename countT>
const bool KMerCountTab<countT>::GetIndexOnDsk(size_t &index_val, const uint64_t kmer_code) const
/* ------------------------------------------------------------------------ *\
    Arg:    k-mer unique code
    Value:  a pair object indicating
                1. size_t for the indexed position on disk (as refArg)
                2. bool for whether the query was succeeded
    Func:   return the k-mer indexed position on disk by k-mer unique code
\* ------------------------------------------------------------------------ */
{
    auto iter = kmer_count_dict_.find(kmer_code);
    if (iter != kmer_count_dict_.cend())
    {
        index_val = iter->second;
        return true;
    }
    else
    {
        return false;
    }
}

#endif //KAMRAT_DATASTRUCT_KMERCOUNTTABTMPLT_HPP