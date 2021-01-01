#ifndef KAMRAT_DATASTRUCT_MERGEKNOT_H
#define KAMRAT_DATASTRUCT_MERGEKNOT_H

#include <string>
#include <map>

class MergeKnot
{
public:
    MergeKnot() noexcept;
    void AddContig(size_t contig_serial, bool is_rc, const std::string &&which_to_set) noexcept;
    const size_t GetSerial(const std::string &&which_to_get) const noexcept;
    const bool IsRC(const std::string &&which_to_get) const noexcept;
    const bool IsMergeable() const noexcept;

private:
    size_t pred_serial_, succ_serial_;
    bool is_pred_rc_, has_pred_, is_succ_rc_, has_succ_, has_ambiguity_;
};

using fix2knot_t = std::map<uint64_t, MergeKnot>; // map for making results not dependent to k-mers' input order, but their sequence order

#endif //KAMRAT_DATASTRUCT_MERGEKNOT_H