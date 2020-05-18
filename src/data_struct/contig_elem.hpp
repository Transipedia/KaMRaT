#ifndef KAMRAT_DATASTRUCT_CONTIGELEM_HPP
#define KAMRAT_DATASTRUCT_CONTIGELEM_HPP

#include <cstdint>

#include "seq_elem.hpp"

class ContigElem : public SeqElem
{
public:
    ContigElem(const std::string &tag,
               const std::string &seq,
               float rep_value,
               uint64_t initial_kmer_code);
    const float GetRepValue() const;
    const std::vector<uint64_t> GetMemberKMerVect() const;
    const size_t GetNbMemberKMer() const;
    const void LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap);
    const void RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap);

private:
    float rep_value_;
    std::vector<uint64_t> member_kmer_vect_;
};

#endif //KAMRAT_DATASTRUCT_CONTIGELEM_HPP