#ifndef KAMRAT_DATASTRUCT_CONTIGLIST_HPP
#define KAMRAT_DATASTRUCT_CONTIGLIST_HPP

#include <string>
#include <cstdint>
#include <vector>

class ContigElem
{
public:
    ContigElem(const std::string &seq, size_t init_serial);

    const bool IsUsed() const;
    const void SetUsed();
    const std::string &GetSeq() const;

    const unsigned int GetNbKMer() const;
    const size_t GetHeadKMerSerial(const bool if_need_reverse) const;
    const size_t GetRepKMerSerial() const;
    const size_t GetRearKMerSerial(const bool if_need_reverse) const;
    const std::vector<size_t> &GetMemKMerSerialVect() const;

    const void LeftExtend(ContigElem &left_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void RightExtend(ContigElem &right_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void ReverseComplement();

private:
    bool is_used_;                             // if the contig is used for merge
    std::string seq_;                          // contig sequence
    size_t rep_kmer_serial_;                   // representative k-mer serial in k-mer count table
    std::vector<size_t> mem_kmer_serial_vect_; // member k-mers for calculating mean count of all member k-mers,
                                               // since seq_ is not necessary for doing this (some k-mer maybe skipped during intervention)
                                               // head k-mer as the first, rear k-mer as the last, others whatever order
};

// using code2contig_t = std::map<uint64_t, ContigElem>;
using contigvect_t = std::vector<ContigElem>;

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP