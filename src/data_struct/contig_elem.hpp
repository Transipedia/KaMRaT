#ifndef KAMRAT_DATASTRUCT_CONTIGLIST_HPP
#define KAMRAT_DATASTRUCT_CONTIGLIST_HPP

#include <string>
#include <cstdint>
#include <set>
#include <map>

#include "seq_elem.hpp"

class ContigElem : public SeqElem
{
public:
    ContigElem(const std::string &seq, float rep_value, size_t init_serial);
    const bool IsUsed() const;
    const void SetUsed();
    const size_t GetHeadSerial() const;
    const size_t GetRearSerial() const;
    const size_t GetRepSerial() const;
    const float GetRepValue() const;
    const std::set<size_t> &GetKMerSerialSet() const;
    const size_t GetNbMemberKMer() const;
    const void LeftMerge(const ContigElem &left_contig_elem, unsigned int n_overlap);
    const void RightMerge(const ContigElem &right_contig_elem, unsigned int n_overlap);
    const void SelfReverseComplement();

private:
    bool is_used_;
    const size_t rep_serial_;
    size_t head_serial_, rear_serial_;
    const float rep_value_;
    std::set<size_t> kmer_serial_set_;
};

using code2contig_t = std::map<uint64_t, ContigElem>;

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP