#ifndef KAMRAT_DATASTRUCT_SEQELEM_HPP
#define KAMRAT_DATASTRUCT_SEQELEM_HPP

#include <string>
#include <map>

class SeqElem
{
public:
    SeqElem(const std::string &seq);
    const std::string GetSeq() const;

protected:
    std::string seq_;
};

using tag2seq_t = std::map<std::string, SeqElem>;

#endif //KAMRAT_DATASTRUCT_SEQELEM_HPP
