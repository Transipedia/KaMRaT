#ifndef KAMRAT_DATASTRUCT_SEQELEM_HPP
#define KAMRAT_DATASTRUCT_SEQELEM_HPP

#include <vector>
#include <unordered_set>
#include <unordered_map>

class SeqElem
{
public:
    SeqElem(const std::string &tag, const std::string &seq);
    const bool IsAccessible() const;
    void MakeUnaccessible();
    const std::string GetTag(const std::string &tag_type) const;
    const std::string GetSeq() const;
    const std::string GetHeadKMer(unsigned int k_len) const;
    const std::string GetRearKMer(unsigned int k_len) const;

private:
    bool is_accessible_; // whether the contig is accessable
    std::string name_;   // sequence name
    std::string seq_;    // sequence
};

using seqVect_t = std::vector<SeqElem>;

#endif //KAMRAT_DATASTRUCT_SEQELEM_HPP
