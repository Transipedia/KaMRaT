#ifndef KAMRAT_MERGE_OVERLAPKNOT_H
#define KAMRAT_MERGE_OVERLAPKNOT_H

#include <unordered_map>

class OverlapKnot
{
public:
    OverlapKnot();
    void SetPredtig(size_t contig_Serial, bool is_rc);
    void SetSucctig(size_t contig_Serial, bool is_rc);
    const size_t GetPredtigSerial() const;
    const size_t GetSucctigSerial() const;
    const bool IsPredtigRC() const;
    const bool IsSucctigRC() const;
    const bool IsSingleKnot() const;
private:
    size_t predtig_serial_, succtig_serial_;
    bool predtig_rc_, succtig_rc_, has_ambiguity_;
};

#endif //KAMRATMERGE_OVERLAPKNOT_H