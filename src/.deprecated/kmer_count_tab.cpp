#include <sstream>
#include <limits> // numeric_limits<streamsize>

#include "kmer_count_tab.hpp"

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded);

KMerCountTab::KMerCountTab(std::vector<std::pair<std::string, size_t>> &kmer_vect, const std::string &idx_info_path)
{
    std::string term;
    std::ifstream idx_info(idx_info_path);
    if (!idx_info.is_open())
    {
        throw std::invalid_argument("KaMRaT index folder not found or may be corrupted");
    }
    
    idx_info >> nb_smp_ >> k_len_ >> term; // load: sample number, k-mer length, and strandedness
    stranded_ = (term == "T");

    colname_vect_.resize(nb_smp_ + 1); // column names include the first column
    for (size_t i(0); i <= nb_smp_; ++i)
    {
        idx_info >> colname_vect_[i]; // load: column name vector
    }

    idx_info.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignore: sum vector position

    for (size_t s(0), pos; idx_info >> term >> pos; ++s) // load: {feature, position} pair
    {
        pos_vect_.emplace_back(pos);
        kmer_serial_map_.insert({Seq2Int(term, k_len_, stranded_), s});
        kmer_vect.emplace_back(std::make_pair(term, pos));
    }

    idx_info.close();
}

const size_t KMerCountTab::GetSampleNumber() const
{
    return nb_smp_;
}

const size_t KMerCountTab::GetKLen() const
{
    return k_len_;
}

const bool KMerCountTab::IsStranded() const
{
    return stranded_;
}

const std::map<uint64_t, size_t> &KMerCountTab::GetKMerSerialMap() const
{
    return kmer_serial_map_;
}

// const std::vector<float> &KMerCountTab::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t serial) const
// {
//     count_vect.resize(nb_smp_);
//     idx_mat.seekg(rowpos_vect_[serial]);
//     idx_mat.read(reinterpret_cast<char *>(&count_vect[0]), nb_smp_ * sizeof(float));
//     return count_vect;
// }

// const std::vector<double> &KMerCountTab::GetSumVect(std::vector<double> &sum_vect, std::ifstream &idx_mat) const
// {
//     sum_vect.resize(nb_smp_);
//     idx_mat.seekg(smp_sum_pos_);
//     idx_mat.read(reinterpret_cast<char *>(&sum_vect[0]), nb_smp_ * sizeof(double)); // load: normalization factor
//     return sum_vect;
// }

// const std::vector<float> &KMerElem::GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_file, const size_t nb_count,
//                                                  const std::vector<double> &nf_vect) const
// {
//     count_vect.resize(nb_count);
//     idx_file.seekg(idx_pos_);
//     idx_file.read(reinterpret_cast<char *>(&count_vect[0]), nb_count * sizeof(float));
//     for (size_t i(0); i < nb_count; ++i)
//     {
//         count_vect[i] *= nf_vect[i];
//     }
//     return count_vect;
// }
