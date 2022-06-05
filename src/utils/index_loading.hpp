

#ifndef ILOAD_H
#define ILOAD_H

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, const std::string &idx_meta_path);

const std::vector<double> &ComputeNF(std::vector<double> &smp_sum_vect, const size_t nb_smp);

void LoadPosVect(std::vector<size_t> &pos_vect, const std::string &idx_pos_path, const bool need_skip_code);

const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat,
                                       const size_t pos, const size_t nb_smp);

#endif