#ifndef KAMRAT_DATASTRUCT_CONTIGLIST_HPP
#define KAMRAT_DATASTRUCT_CONTIGLIST_HPP

#include <cstdint>
#include <set>
#include <map>

#include "seq_elem.hpp"

class ContigElem : public SeqElem
{
public:
    ContigElem(const std::string &seq, float rep_value, uint64_t initial_kmer_code);
    const float GetRepValue() const;
    const std::set<uint64_t> GetMemberKMerSet() const;
    const size_t GetNbMemberKMer() const;
    const bool LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap);
    const bool RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap);

private:
    float rep_value_;
    std::set<uint64_t> member_kmer_set_;
};

using code2contig_t = std::map<uint64_t, ContigElem>;

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP