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
    ContigElem(const std::string &seq, float rep_value, uint64_t init_uniqcode, size_t init_serial);
    const bool IsUsed() const;
    const void SetUsed();
    const unsigned int GetNbKMer() const;
    const size_t GetHeadKMerSerial(const bool if_need_reverse) const;
    const size_t GetRearKMerSerial(const bool if_need_reverse) const;
    const void LeftExtend(const ContigElem &left_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void RightExtend(const ContigElem &right_contig_elem, bool need_right_rc, unsigned int n_overlap);

private:
    bool is_used_;
    unsigned int nb_kmer_;
    size_t head_kmer_serial_, rear_kmer_serial_;
};

using code2contig_t = std::map<uint64_t, ContigElem>;
using contigvect_t = std::vector<ContigElem>;

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP