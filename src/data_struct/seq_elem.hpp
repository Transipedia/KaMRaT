#ifndef KAMRAT_DATASTRUCT_SEQELEM_HPP
#define KAMRAT_DATASTRUCT_SEQELEM_HPP

#include <cstdint>
#include <vector>
#include <unordered_set>
#include <unordered_map>

class SeqElem
{
public:
    SeqElem(const std::string &tag, const std::string &seq);
    const bool IsAccessible() const;
    void MakeUnaccessible();
    const std::string GetTag() const;
    const std::string GetSeq() const;
    const std::string GetKMerAt(size_t pos, unsigned int k_len) const;
    const std::string GetHeadKMer(unsigned int k_len) const;
    const std::string GetRearKMer(unsigned int k_len) const;

private:
    bool is_accessible_; // whether the contig is accessable
    std::string tag_;    // sequence name
    std::string seq_;    // sequence
};

using seqVect_t = std::vector<SeqElem>;

#endif //KAMRAT_DATASTRUCT_SEQELEM_HPP