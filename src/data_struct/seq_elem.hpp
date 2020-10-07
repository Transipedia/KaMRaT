#ifndef KAMRAT_DATASTRUCT_SEQELEM_HPP
#define KAMRAT_DATASTRUCT_SEQELEM_HPP

#include <string>
#include <map>
#include <vector>
#include <cstdint>

class SeqElem
{
public:
    SeqElem(const std::string &seq, uint64_t uniq_code, float score);
    const std::string GetSeq() const;
    const uint64_t GetUniqCode() const;
    const float GetScore(const std::string &&mode) const;
    void ScaleScore(float fact, float lower_lim, float upper_lim);

protected:
    uint64_t uniq_code_; // could be any unique code, either serial number or hashed representative k-mer code
    std::string seq_;    // sequence
    float score_;        // score related to the sequence, unmodifiable
    float final_score_;  // a modifiable score for p-value adjustement
};

using name2seq_t = std::map<std::string, SeqElem>;
using seqVect_t = std::vector<SeqElem>;

#endif //KAMRAT_DATASTRUCT_SEQELEM_HPP
