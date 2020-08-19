#ifndef KAMRAT_MERGE_CONTIGINFO_H
#define KAMRAT_MERGE_CONTIGINFO_H

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>

class ContigInfo
{
public:
    ContigInfo(const std::string &kmer_seq, double rep_value, size_t count_pos);
    ContigInfo(const std::string &kmer_seq, double rep_value, const std::vector<double> &counts);
    ContigInfo(const ContigInfo &predtig, const ContigInfo &succtig, size_t n_overlap, const std::string &quant_mode);
    const bool IsUsed() const;
    void SetUsed();
    const std::string &GetTag() const;
    const std::string &GetSequence() const;
    void ReverseComplement();
    const size_t GetNbKmer() const;
    const double GetRepValue() const;
    const std::vector<size_t> &GetCountPos() const;
    const std::vector<double> &GetCounts() const;
    const std::vector<double> &GetHeadCounts() const;
    const std::vector<double> &GetRearCounts() const;
    const size_t GetHeadPos() const;
    const size_t GetRearPos() const;

private:
    bool is_used_;
    std::string tag_;               // name for kamratReduce or tag for kamratMerge
    std::string seq_;               // feature sequence
    size_t nb_kmer_;                // member k-mer number
    double rep_value_;              // the lower the better
    std::vector<size_t> count_pos_; // source for counting,
                                    // contains rep-kmer code in rep-kmer mode,
                                    // contains all member k-mers in mean mode
    size_t head_pos_, rear_pos_;    // count position for the head and rear k-mers
    std::vector<double> counts_, head_counts_, rear_counts_; // for in-memory mode, to deprecate
};

#endif //KAMRAT_MERGE_CONTIGINFO_H
