#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <limits>
#include <ctime>

#include "runinfo_files/filter_runinfo.hpp"

#define RESET "\033[0m"
#define BOLDYELLOW "\033[1m\033[33m"

void LoadIndexMeta(size_t &nb_smp_all, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, const std::string &idx_meta_path);                // in utils/index_loading.cpp
const std::vector<double> &ComputeNF(std::vector<double> &smp_sum_vect, const size_t nb_smp);                // in utils/index_loading.cpp
void LoadPosVect(std::vector<size_t> &pos_vect, const std::string &idx_pos_path, const bool need_skip_code); // in utils/index_loading.cpp
const std::vector<float> &GetCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat,
                                       const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp

const std::pair<size_t, size_t> ParseDesign(std::vector<bool> &filter_stat_vect, const std::string &dsgn_path,
                                            const std::vector<std::string> &colname_vect, const size_t nb_smp)
{
    std::ifstream dsgn_file(dsgn_path);
    if (!dsgn_file.is_open())
    {
        throw std::invalid_argument("error open design file: " + dsgn_path);
    }
    std::string smp_name, filter_stat;
    std::unordered_map<std::string, bool> smp_isup_map; // false means DOWN, true means UP
    while (std::getline(dsgn_file, smp_name))
    {
        size_t split_pos = smp_name.find_first_of("\t ");
        if (split_pos != std::string::npos)
        {
            filter_stat = smp_name.substr(split_pos + 1);
            smp_name = smp_name.substr(0, split_pos);
        }
        if (filter_stat == "UP")
        {
            smp_isup_map.insert({smp_name, true});
        }
        else if (filter_stat == "DOWN")
        {
            smp_isup_map.insert({smp_name, false});
        }
        else
        {
            throw std::invalid_argument("unknown filter status " + filter_stat + ", only UP or DOWN are expected");
        }
    }
    filter_stat_vect.resize(nb_smp);
    size_t nb_smp_up(0), nb_smp_down(0);
    for (size_t i(1); i <= nb_smp; ++i)
    {
        const auto &it = smp_isup_map.find(colname_vect[i]);
        if (it != smp_isup_map.cend())
        {
            filter_stat_vect[i - 1] = it->second;
            it->second ? nb_smp_up++ : nb_smp_down++;
        }
    }
    // for (size_t i(0); i < nb_smp; ++i)
    // {
    //     std::cout << colname_vect[i + 1] << "\t" << (filter_stat_vect[i] ? "UP" : "DOWN") << std::endl;
    // }
    dsgn_file.close();
    return std::make_pair(nb_smp_up, nb_smp_down);
}

void ScanPrint(std::ifstream &idx_mat, const std::vector<size_t> &ft_pos_vect, const std::vector<bool> &filter_stat_vect,
               const size_t up_min_abd, const size_t up_min_rec, const size_t down_max_abd, const size_t down_min_rec,
               const size_t nb_smp, const bool reverse_filter, const bool with_counts)
{
    std::vector<float> count_vect;
    std::string ft_name;
    for (size_t ft_pos : ft_pos_vect)
    {
        GetCountVect(count_vect, idx_mat, ft_pos, nb_smp);
        idx_mat >> ft_name;
        size_t up_rec(0), down_rec(0);
        for (size_t i_smp(0); i_smp < nb_smp; ++i_smp)
        {
            if (filter_stat_vect[i_smp] && count_vect[i_smp] >= up_min_abd)
            {
                ++up_rec;
            }
            else if (!filter_stat_vect[i_smp] && count_vect[i_smp] <= down_max_abd)
            {
                ++down_rec;
            }
        }
        if (reverse_filter != (up_rec >= up_min_rec && down_rec >= down_min_rec)) // !reverse && eligible || reverse && !eligible
        {
            std::cout << ft_name;
            if (with_counts)
            {
                for (const float x : count_vect)
                {
                    std::cout << "\t" << x;
                }
            }
            else
            {
                std::cout << "\t0\t1\t";
                std::cout.write(reinterpret_cast<char *>(&ft_pos), sizeof(size_t));
            }
            std::cout << std::endl;
        }
    }
}

int FilterMain(int argc, char *argv[])
{
    FilterWelcome();

    std::clock_t begin_time = clock();
    std::string idx_dir, dsgn_path, out_path;
    size_t up_min_rec(0), up_min_abd(0), down_min_rec(0), down_max_abd(std::numeric_limits<size_t>::max()), nb_smp, k_len;
    bool reverse_filter(false), with_counts(false), _stranded; // _stranded not needed
    std::vector<std::string> colname_vect;

    ParseOptions(argc, argv, idx_dir, dsgn_path, up_min_abd, up_min_rec, down_max_abd, down_min_rec, reverse_filter, out_path, with_counts);
    PrintRunInfo(idx_dir, dsgn_path, up_min_abd, up_min_rec, down_max_abd, down_min_rec, reverse_filter, out_path, with_counts);
    LoadIndexMeta(nb_smp, k_len, _stranded, colname_vect, idx_dir + "/idx-meta.bin");

    std::vector<bool> filter_stat_vect;
    const std::pair<size_t, size_t> &&dsgn_info = ParseDesign(filter_stat_vect, dsgn_path, colname_vect, nb_smp);

    if (dsgn_info.first < up_min_rec)
    {
        std::cerr << BOLDYELLOW << "[warning] " << RESET << "UP column number smaller than given minimum recurrence threshold: "
                  << dsgn_info.first << "<" << up_min_rec << std::endl
                  << std::endl;
    }
    if (dsgn_info.second < down_min_rec)
    {
        std::cerr << BOLDYELLOW << "[warning] " << RESET << "DOWN column number smaller than given minimum recurrence threshold: "
                  << dsgn_info.second << "<" << down_min_rec << std::endl
                  << std::endl;
    }
    if (with_counts)
    {
        std::cout << colname_vect[0];
        for (size_t i_col(1); i_col <= nb_smp; ++i_col)
        {
            std::cout << "\t" << colname_vect[i_col];
        }
        std::cout << std::endl;
    }
    std::vector<size_t> ft_pos_vect;
    LoadPosVect(ft_pos_vect, idx_dir + "/idx-pos.bin", k_len != 0);
    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    std::ofstream out_file;
    if (!out_path.empty())
    {
        out_file.open(out_path);
        if (!out_file.is_open())
        {
            throw std::domain_error("cannot open file: " + out_path);
        }
    }
    auto backup_buf = std::cout.rdbuf();
    if (!out_path.empty()) // output to file if a path is given, to screen if not
    {
        std::cout.rdbuf(out_file.rdbuf());
    }
    ScanPrint(idx_mat, ft_pos_vect, filter_stat_vect, up_min_abd, up_min_rec, down_max_abd, down_min_rec, nb_smp, reverse_filter, with_counts);
    idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}
