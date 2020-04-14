#ifndef KAMRAT_QUERY_CONTIGSAMPLECOUNT_HPP
#define KAMRAT_QUERY_CONTIGSAMPLECOUNT_HPP

#include <string>
#include <vector>

class ContigSampleCount
{
public:
    void AddKmerCounts(const std::vector<size_t> &sample_count_one_kmer);

private:
    std::vector<std::vector<size_t>> sample_count_all_kmers_;
};

#endif //KAMRAT_QUERY_CONTIGSAMPLECOUNT_HPP