#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &seq)
    : seq_(seq)
{
}

const std::string SeqElem::GetSeq() const
{
    return seq_;
}
