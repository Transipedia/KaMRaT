#include "seq_elem.hpp"

SeqElem::SeqElem(const std::string &name, const std::string &seq)
    : is_accessible_(true),
      name_(name),
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

const std::string SeqElem::GetTag(const std::string &tag_type) const
{
    if (tag_type == "name")
    {
        return name_;
    }
    else if(tag_type == "seq")
    {
        return GetSeq();
    }
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