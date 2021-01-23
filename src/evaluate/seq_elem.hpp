#ifndef KAMRAT_EVALUATE_SEQELEM_HPP
#define KAMRAT_EVALUATE_SEQELEM_HPP

#include <string>
#include <vector>
#include <memory>

class SeqElem
{
public:
    SeqElem(const std::string &name, const std::string &seq);

    const std::string &GetName() const;
    const std::string &GetSeq() const;

private:
    const std::string name_, seq_; // sequence name and sequence
};

using fastaVect_t = std::vector<std::unique_ptr<SeqElem>>;

#endif //KAMRAT_EVALUATE_SEQELEM_HPP