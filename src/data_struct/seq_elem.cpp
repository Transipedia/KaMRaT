#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &tag, const std::string &seq)
    : is_accessible_(true),
      tag_(tag),
      seq_(seq)
{
}

const bool SeqElem::IsAccessible() const
{
    return is_accessible_;
}

void SeqElem::MakeUnaccessible()
{
    is_accessible_ = false;
}

const std::string SeqElem::GetTag() const
{
    return tag_;
}

const std::string SeqElem::GetSeq() const
{
    return seq_;
}

