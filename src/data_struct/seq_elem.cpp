#include <algorithm>

#include "seq_elem.hpp"

inline void ToComplement(std::string &seq)
{
    auto lambda_trans = [](const char c) {
        switch (c)
        {
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        default:
            return 'N';
        }
    };
    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda_trans);
}

SeqElem::SeqElem(const std::string &seq)
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

const void SeqElem::SelfReverseComplement()
{
    reverse(seq_.begin(), seq_.end());
    ToComplement(seq_);
}