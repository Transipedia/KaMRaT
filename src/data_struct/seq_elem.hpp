#ifndef KAMRAT_DATASTRUCT_SEQELEM_HPP
#define KAMRAT_DATASTRUCT_SEQELEM_HPP

#include <string>
#include <map>
#include <vector>

class SeqElem
{
public:
    SeqElem(const std::string &seq, const size_t serial, const float score);
    const std::string GetSeq() const;
    const size_t GetSerial() const;
    const float GetScore(const std::string &&mode) const;
    void ScaleScore(float fact, float lower_lim, float upper_lim);

protected:
    std::string seq_;
    const size_t serial_;
    const float score_;
    float final_score_; // a modifiable score for p-value adjustement
};

using tag2seq_t = std::map<std::string, SeqElem>;
using seqVect_t = std::vector<SeqElem>;

#endif //KAMRAT_DATASTRUCT_SEQELEM_HPP
