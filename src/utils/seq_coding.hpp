#ifndef KAMRAT_UTILS_SEQCODING_HPP
#define KAMRAT_UTILS_SEQCODING_HPP

#include <cstdint>
#include <cstddef>
#include <string>
#include <stdexcept>

static inline uint8_t Nuc2Num(const char nuc)
{
    switch (nuc)
    {
    case 'a':
    case 'A':
        return 0;
    case 'c':
    case 'C':
        return 1;
    case 'g':
    case 'G':
        return 2;
    case 't':
    case 'T':
        return 3;
    default:
        return 255;
    }
}

static inline char Num2Nuc(const uint8_t c)
{
    switch (c)
    {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        return 'N';
    }
}

static inline uint64_t GetRC(const uint64_t code, size_t k_length)
{
    uint64_t mask, tmp_code(code);
    if (k_length == 32)
        mask = ~0;
    else
        mask = ((uint64_t)1 << (2 * k_length)) - 1;

    tmp_code ^= mask;

    uint64_t mask_lsb;
    mask_lsb = 3; // corresponds to the rightmost nucleotide //
    uint64_t shift = 0;
    uint64_t result = 0;
    for (size_t j = 0; j < k_length; ++j)
    {
        result <<= 2;
        result |= (tmp_code & mask_lsb) >> shift; // get the leftmost nucleotide and put it at the end //
        mask_lsb <<= 2;
        shift += 2;
    }

    return result;
}

static inline void Int2Seq(std::string &seq, const uint64_t code, const size_t k_length)
{
    uint64_t tmp_code(code);
    seq.resize(k_length);
    uint64_t mask = 3;
    for (size_t i = 0; i < k_length; ++i)
    {
        seq[k_length - i - 1] = Num2Nuc(tmp_code & mask);
        tmp_code >>= 2;
    }
}

static inline uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded)
{
    uint64_t code = 0;
    for (size_t i = 0; i < k_length; ++i)
    {
        code <<= 2;
        code |= Nuc2Num(seq[i]);
    }
    if (!stranded)
    {
        uint64_t code_rc(GetRC(code, k_length));
        if (code_rc < code)
        {
            code = code_rc;
        }
    }
    return code;
}

static inline uint64_t MutIntAtPos(uint64_t code, const size_t k_length, const int pos, const char n)
{
    uint64_t x = (k_length - pos - 1) * 2;
    uint64_t up = 1;
    if (n == 'A')
    {
        code &= ~(up << x);
        code &= ~(up << (x + 1));
    }
    else if (n == 'C')
    {
        code |= up << x;
        code &= ~(up << (x + 1));
    }
    else if (n == 'G')
    {
        code &= ~(up << x);
        code |= up << (x + 1);
    }
    else
    {
        code |= up << x;
        code |= up << (x + 1);
    }
    return code;
}

static inline uint64_t NextSeq(uint64_t code, const size_t k_length, const char new_nuc, const bool stranded)
{
    uint64_t x = (k_length - 1) * 2, up = 3;
    // Unset the first nucleotide
    code &= ~(up << x);
    code <<= 2;
    code = MutIntAtPos(code, k_length, k_length - 1, new_nuc);
    if (!stranded)
    {
        uint64_t code_rc = GetRC(code, k_length);
        if (code_rc < code)
        {
            code = code_rc;
        }
    }
    return code;
}

#endif //KAMRAT_UTILS_SEQCODING_HPP