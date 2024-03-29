#include <string>
#include <limits>

inline uint8_t Nuc2Num(const char nuc)
{
    if (nuc == 'a' || nuc == 'A')
    {
        return 0;
    }
    else if (nuc == 'c' || nuc == 'C')
    {
        return 1;
    }
    else if (nuc == 'g' || nuc == 'G')
    {
        return 2;
    }
    else if (nuc == 't' || nuc == 'T')
    {
        return 3;
    }
    else
    {
        return std::numeric_limits<uint8_t>::max();
    }
}

inline char Num2Nuc(const uint8_t c)
{
    if (c == 0)
    {
        return 'A';
    }
    else if (c == 1)
    {
        return 'C';
    }
    else if (c == 2)
    {
        return 'G';
    }
    else if (c == 3)
    {
        return 'T';
    }
    else
    {
        return 'N';
    }
}

uint64_t GetRC(const uint64_t code, size_t k_length)
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

void Int2Seq(std::string &seq, const uint64_t code, const size_t k_length)
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

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded)
{
    uint64_t code = 0;
    for (size_t i = 0; i < k_length; ++i)
    {
        code <<= 2;
        code |= Nuc2Num(seq[i]);
    }
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

inline uint64_t MutIntAtPos(uint64_t code, const size_t k_length, const int pos, const char n)
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

uint64_t NextCode(uint64_t code, const size_t k_length, const char new_nuc)
// Attention: should not return the RC code if unstranded !
//            returning of RC code would disturb the calculation of downstreaming k-mers.
{
    uint64_t x = (k_length - 1) * 2, up = 3;
    // Unset the first nucleotide
    code &= ~(up << x); // TO TEST: is this necessary ?
    code <<= 2;
    code = MutIntAtPos(code, k_length, k_length - 1, new_nuc);
    return code;
}
