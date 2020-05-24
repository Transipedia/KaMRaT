#ifndef KAMRAT_DATASTRUCT_MERGEKNOT_H
#define KAMRAT_DATASTRUCT_MERGEKNOT_H

#include <cstdint>
#include <string>
#include <map>

class MergeKnot
{
public:
    MergeKnot();
    void AddContig(uint64_t kmer_code, bool is_rc, const std::string &which_to_set); // representative k-mer code for contig calling
    const uint64_t GetPredCode() const;
    const uint64_t GetSuccCode() const;
    const bool IsPredRC() const;
    const bool IsSuccRC() const;
    const bool IsMergeable() const;
private:
    uint64_t pred_code_, succ_code_;
    bool is_pred_rc_, is_succ_rc_, has_pred_, has_succ_, has_ambiguity_;
};

using fix2knot_t = std::map<uint64_t, MergeKnot>;

#endif //KAMRAT_DATASTRUCT_MERGEKNOT_H