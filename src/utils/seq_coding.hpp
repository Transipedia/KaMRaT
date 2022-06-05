
#ifndef SEQ_H
#define SEQ_H

const uint64_t Seq2Int(const std::string &seq, const size_t k_len, const bool stranded);
uint64_t GetRC(const uint64_t code, size_t k_length);
uint64_t NextCode(uint64_t code, const size_t k_length, const char new_nuc);


#endif
