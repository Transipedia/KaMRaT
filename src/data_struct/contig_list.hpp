#ifndef KAMRAT_DATASTRUCT_CONTIGLIST_HPP
#define KAMRAT_DATASTRUCT_CONTIGLIST_HPP

#include <cstdint>
#include <vector>

#include "seq_elem.hpp"

class ContigList
{
public:

private:
    seqVect_t seq_list_;
    
};


// {
// public:
//     ContigElem(const std::string &tag,
//                const std::string &seq,
//                float rep_value,
//                uint64_t initial_kmer_code);
//     const float GetRepValue() const;
//     const std::vector<uint64_t> GetMemberKMerVect() const;
//     const size_t GetNbMemberKMer() const;
//     const bool LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap);
//     const bool RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap);

// private:
//     float rep_value_;
//     std::vector<uint64_t> member_kmer_vect_;
// };

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP