#ifndef KAMRAT_UTILS_DNA_SEQUENCE_HPP
#define KAMRAT_UTILS_DNA_SEQUENCE_HPP

#include <string>
#include <stdexcept>
#include <algorithm>

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
            throw std::domain_error("Invalid nucleotide: " + c);
        }
    };
    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda_trans);
}

inline void ToReverseComplement(std::string &seq)
{
    reverse(seq.begin(), seq.end());
    ToComplement(seq);
}

#endif //KAMRAT_UTILS_DNA_SEQUENCE_HPP