
#ifndef SEQ_H
#define SEQ_H

const uint64_t Seq2Int(const std::string &seq, const size_t k_len, const bool stranded);
uint64_t GetRC(const uint64_t code, size_t k_length);

#endif
