#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &name, const std::string &seq)
    : seq_(seq)
{
}

const std::string SeqElem::GetSeq() const
{
    return seq_;
}

const std::string SeqElem::GetHeadKMer(const unsigned int k_len) const
{
    return seq_.substr(0, k_len);
}

const std::string SeqElem::GetRearKMer(const unsigned int k_len) const
{
    return seq_.substr(seq_.size() - k_len);
}