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
    const std::vector<size_t> &GetMemKMerSerialVect(std::vector<size_t> &mem_kmer_vect) const;
    
    const void LeftExtend(const ContigElem &left_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void RightExtend(const ContigElem &right_contig_elem, bool need_right_rc, unsigned int n_overlap);

private:
    bool is_used_;
    std::string seq_;
    size_t head_kmer_serial_, rep_kmer_serial_, rear_kmer_serial_;
    std::vector<size_t> mem_kmer_serial_vect_; // for storing contig's member k-mers, used for calculating mean count of contigs (sequence is not sufficient)
};

// using code2contig_t = std::map<uint64_t, ContigElem>;
using contigvect_t = std::vector<ContigElem>;

#endif //KAMRAT_DATASTRUCT_CONTIGLIST_HPP