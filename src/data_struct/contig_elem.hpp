#ifndef KAMRAT_MERGE_CONTIGELEM_HPP
#define KAMRAT_MERGE_CONTIGELEM_HPP

#include <string>
#include <vector>
#include <memory>
#include <fstream>

class ContigElem
{
public:
    ContigElem(const std::string &seq, size_t pos, float val);

    const std::string &GetSeq() const;
    const size_t GetRepPos() const;
    const float GetRepVal() const;
    const size_t GetHeadPos(bool need_reverse) const;
    const size_t GetRearPos(bool need_reverse) const;
    const size_t GetNbMemKmer() const;
    const std::vector<size_t> &GetMemPosVect() const;
    const void LeftExtend(std::unique_ptr<ContigElem> left_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void RightExtend(std::unique_ptr<ContigElem> right_contig_elem, bool need_right_rc, unsigned int n_overlap);
    const void ReverseComplement();

private:
    std::string seq_;                  // contig sequence
    const size_t rep_pos_;             // representative pos
    const float rep_val_;              // value for indicating representativeness
    std::vector<size_t> mem_pos_vect_; // component k-mer positions in index
                                       // head pos at first and rear pos at last, others in middle (rep pos repeated)
};

#endif //KAMRAT_MERGE_CONTIGELEM_HPP