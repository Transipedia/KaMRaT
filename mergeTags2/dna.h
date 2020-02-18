#ifndef DNA_H
#define DNA_H

#include <stdint.h>

static const int base_to_int[255] =
    {
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*   0 -   9 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  10 -  19 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  20 -  29 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  30 -  39 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  40 -  49 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  50 -  59 */
        -1, -1, -1, -1, -1, 0, -1, 1, -1, -1,   /*  60 -  69 */
        -1, 2, -1, -1, -1, -1, -1, -1, -1, -1,  /*  70 -  79 */
        -1, -1, -1, -1, 3, 3, -1, -1, -1, -1,   /*  80 -  89 */
        -1, -1, -1, -1, -1, -1, -1, 0, -1, 1,   /*  90 -  99 */
        -1, -1, -1, 2, -1, -1, -1, -1, -1, -1,  /* 100 - 109 */
        -1, -1, -1, -1, -1, -1, 3, 3, -1, -1,   /* 110 - 119 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 120 - 129 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 130 - 139 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 140 - 149 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 150 - 159 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 160 - 169 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 170 - 179 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 180 - 189 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 190 - 199 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 200 - 209 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 210 - 219 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 220 - 229 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 230 - 239 */
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 240 - 249 */
        -1, -1, -1, -1, -1                      /* 250 - 254 */
};

static const char NUCLEOTIDES[4] = {'A', 'C', 'G', 'T'};

static const int NB_NUCLEOTIDES = 4;

static inline void int_to_dna(uint64_t code, size_t dna_length, char *dna)
{
    uint64_t mask = 3;
    for (size_t i = 0; i < dna_length; i++)
    {
        dna[dna_length - i - 1] = NUCLEOTIDES[code & mask];
        code >>= 2;
    }
}

static inline uint64_t dna_to_int(const char *dna, size_t dna_length)
{
    // TODO we should check that dna_length is smaller that 32.
    uint64_t code = 0;
    for (size_t i = 0; i < dna_length; i++)
    {
        code <<= 2;
        code |= base_to_int[(int)dna[i]];
        //fprintf(stderr, "%c => %d => %" PRIu64 "\n", dna[i], base_to_int[(int)dna[i]], code);
    }
    return code;
}

static inline uint64_t int_revcomp(uint64_t factor, uint32_t length)
{
    uint64_t mask;
    if (length == 32)
        mask = ~0;
    else
        mask = ((uint64_t)1 << (2 * length)) - 1;

    factor ^= mask;

    uint64_t mask_lsb;
    // Corresponds to the rightmost nucleotide
    mask_lsb = 3;
    uint64_t shift = 0;
    uint64_t result = 0;
    for (size_t j = 0; j < length; j++)
    {
        result <<= 2;
        // get the leftmost nucleotide and put it at the end
        result |= (factor & mask_lsb) >> shift;
        mask_lsb <<= 2;
        shift += 2;
    }

    return result;
}

static inline uint64_t mut_int_dna(uint64_t code, size_t dna_length, int pos, char n)
{
    uint64_t x = (dna_length - pos - 1) * 2;
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

static inline uint64_t next_kmer(size_t dna_length, char new_nuc, uint64_t code)
{
    uint64_t x = (dna_length - 1) * 2, up = 3;
    // Unset the first nucleotide
    code &= ~(up << x);
    code <<= 2;
    return mut_int_dna(code, dna_length, dna_length - 1, new_nuc);
}

uint64_t hash_kmer(const char *seq, size_t k)
{
    uint64_t return_value = 0;
    while (k--)
    {
        char c = *seq++;
        c &= 6;
        c ^= c >> 1;
        uint64_t val = c >> 1;
        return_value <<= 2;
        return_value |= val;
    }
    return return_value;
}

char *revcomp(const char *forward, char *reverse, size_t len)
{
    for (size_t k = 0; k < len; k++)
    {
        char c = forward[k];
        char magic = c & 2 ? 4 : 21;
        reverse[len - k - 1] = c ^ magic;
    }
    reverse[len] = '\0';
    return reverse;
}

#endif
