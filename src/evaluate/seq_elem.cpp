#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &name, const std::string &seq)
    : name_(name), seq_(seq)
{
}

const std::string &SeqElem::GetName() const
{
    return name_;
}

const std::string &SeqElem::GetSeq() const
{
    return seq_;
}
