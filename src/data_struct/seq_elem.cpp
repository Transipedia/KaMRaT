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

const std::string SeqElem::GetKMerAt(const size_t pos, const unsigned int k_len) const
{
    return seq_.substr(pos, k_len);
}

const std::string SeqElem::GetHeadKMer(const unsigned int k_len) const
{
    return GetKMerAt(0, k_len);
}

const std::string SeqElem::GetRearKMer(const unsigned int k_len) const
{
    return GetKMerAt(seq_.size() - k_len, k_len);
}
