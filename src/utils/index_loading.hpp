#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#ifndef ILOAD_H
#define ILOAD_H

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, const std::string &idx_meta_path);
void LoadPosVect(std::vector<size_t> &pos_vect, const std::string &idx_pos_path, const bool need_skip_code);
void LoadCodePosMap(std::map<uint64_t, size_t> &code_set, const std::string &idx_pos_path);
void LoadFeaturePosMap(std::unordered_map<std::string, size_t> &ft_pos_map, std::ifstream &idx_mat, const std::string &idx_pos_path, const bool need_skip_code, const size_t nb_smp);

const std::vector<double> &ComputeNF(std::vector<double> &smp_sum_vect, const size_t nb_smp);

const std::string &GetTagSeq(std::string &tag_str, std::ifstream &idx_mat, const size_t pos, const size_t nb_smp);
const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t pos, const size_t nb_smp);
const std::vector<float> &GetMeanCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp, const std::vector<size_t> &mem_pos_vect);
const std::vector<float> &GetMedianCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp, const std::vector<size_t> &mem_pos_vect);
#endif
