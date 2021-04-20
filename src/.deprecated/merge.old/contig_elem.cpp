#include <algorithm>

#include "contig_elem.hpp"

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

static inline const void ReverseComplementSeq(std::string &seq)
{
    reverse(seq.begin(), seq.end());
    ToComplement(seq);
}

ContigElem::ContigElem(const std::string &seq, const size_t init_serial)
    : is_used_(false), seq_(seq), rep_kmer_serial_(init_serial)
{
    mem_kmer_serial_vect_.push_back(init_serial);
}

const bool ContigElem::IsUsed() const
{
    return is_used_;
}

const void ContigElem::SetUsed()
{
    is_used_ = true;
}

const std::string &ContigElem::GetSeq() const
{
    return seq_;
}

const unsigned int ContigElem::GetNbKMer() const
{
    return mem_kmer_serial_vect_.size();
}

const size_t ContigElem::GetHeadKMerSerial(const bool if_need_reverse) const
{
    return (if_need_reverse ? mem_kmer_serial_vect_.back() : mem_kmer_serial_vect_.front());
}

const uint64_t ContigElem::GetRepKMerSerial() const
{
    return rep_kmer_serial_;
}

const size_t ContigElem::GetRearKMerSerial(const bool if_need_reverse) const
{
    return (if_need_reverse ? mem_kmer_serial_vect_.front() : mem_kmer_serial_vect_.back());
}

const std::vector<size_t> &ContigElem::GetMemKMerSerialVect() const
{
    return mem_kmer_serial_vect_;
}

const void ContigElem::LeftExtend(ContigElem &left_contig_elem, const bool need_left_rc, unsigned int n_overlap)
{
    const size_t left_nb_mem = left_contig_elem.GetNbKMer(), this_nb_mem = GetNbKMer();
    if (need_left_rc)
    {
        left_contig_elem.ReverseComplement();
    }
    seq_ = left_contig_elem.GetSeq() + seq_.substr(n_overlap);
    mem_kmer_serial_vect_.insert(mem_kmer_serial_vect_.end(),
                                 std::make_move_iterator(left_contig_elem.GetMemKMerSerialVect().begin()),
                                 std::make_move_iterator(left_contig_elem.GetMemKMerSerialVect().end()));
    if (this_nb_mem == 1 && left_nb_mem == 1)
    {
        std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_.back());
    }
    else if (this_nb_mem == 1)
    {
        std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_.back());
        std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_[1]);
    }
    else if (left_nb_mem == 1)
    {
        std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_.back());
        std::swap(mem_kmer_serial_vect_.back(), mem_kmer_serial_vect_[this_nb_mem - 1]);
    }
    else
    {
        std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_[this_nb_mem]);    // change contig's head k-mer as the head k-mer of left_contig_elem
        std::swap(mem_kmer_serial_vect_[this_nb_mem - 1], mem_kmer_serial_vect_.back()); // change contig's rear k-mer as the rear k-mer of this_contig
    }
    left_contig_elem.SetUsed();
}

const void ContigElem::RightExtend(ContigElem &right_contig_elem, const bool need_right_rc, unsigned int n_overlap)
{
    if (need_right_rc)
    {
        right_contig_elem.ReverseComplement();
    }
    seq_ = seq_ + right_contig_elem.GetSeq().substr(n_overlap);
    mem_kmer_serial_vect_.insert(mem_kmer_serial_vect_.end(),
                                 std::make_move_iterator(right_contig_elem.GetMemKMerSerialVect().begin()),
                                 std::make_move_iterator(right_contig_elem.GetMemKMerSerialVect().end()));
    right_contig_elem.SetUsed();
}

const void ContigElem::ReverseComplement()
{
    ReverseComplementSeq(seq_);
    std::swap(mem_kmer_serial_vect_.front(), mem_kmer_serial_vect_.back());
}
