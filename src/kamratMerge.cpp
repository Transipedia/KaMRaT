#include <iostream>
#include <ctime>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm> // std::remove_if

#include "merge/merge_runinfo.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"

#define RESET "\033[0m"
#define BOLDYELLOW "\033[1m\033[33m"

using contigVect_t = std::vector<std::unique_ptr<ContigElem>>;

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp
uint64_t GetRC(const uint64_t code, size_t k_length);                                 // in utils/seq_coding.cpp

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y);  // in utils/vect_opera.cpp
const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y); // in utils/vect_opera.cpp
const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y);      // in utils/vect_opera.cpp

void LoadIndexMeta(size_t &nb_smp, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path); // in utils/index_loading.cpp
void LoadCodePosValMap(std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                       std::unordered_map<uint64_t, float> &sel_code_val_map,
                       const std::string &idx_pos_path);                                                              // in utils/index_loading.cpp
const std::string &GetTagSeq(std::string &tag_str, std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp
const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, size_t pos, size_t nb_smp); // in utils/index_loading.cpp
const std::vector<float> &GetMeanCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                           const std::string &seq, const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                                           const size_t k_len, const bool stranded);
const std::vector<float> &GetMedianCountVect(std::vector<float> &count_vect, std::ifstream &idx_mat, const size_t nb_smp,
                                             const std::string &seq, const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                                             const size_t k_len, const bool stranded);

// ----- when the -select argument is provided: merging on subset of index ----- //
const bool LoadSelect(std::unordered_map<uint64_t, float> &sel_kmer_map, const std::string &sel_path, const std::string &rep_mode,
                      const size_t k_len, const bool stranded)
{
    std::ifstream sel_file(sel_path);
    if (!sel_file.is_open())
    {
        throw std::invalid_argument("path for input sequences not accessable: " + sel_path);
    }
    std::string line;
    size_t split_pos;
    float rep_val(0);
    bool has_value(false);
    while (std::getline(sel_file, line))
    {
        split_pos = line.find_first_of("\t ");
        if (split_pos != std::string::npos)
        {
            has_value = true;
            rep_val = std::stof(line.substr(split_pos + 1));
            line = line.substr(0, split_pos);
        }
        if (!sel_kmer_map.insert({Seq2Int(line, k_len, stranded), rep_val}).second) // check unicity of given k-mers
        {
            throw std::domain_error("checking unicity failed: k-mer not unique in the given sequence list: " + line);
        }
    }
    sel_file.close();
    // for (const auto &elem : sel_kmer_map)
    // {
    //     std::cout << elem.first << "\t" << elem.second << std::endl;
    // }
    return has_value;
}

void InitializeContigList(contigVect_t &ctg_vect, std::ifstream &idx_mat, const size_t nb_smp,
                          const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map)
{
    std::string seq;
    for (const auto &elem : code_posval_map)
    {
        idx_mat.seekg(elem.second.first + nb_smp * sizeof(float)); // skip the indexed count vector
        idx_mat >> seq;
        ctg_vect.emplace_back(std::make_unique<ContigElem>(seq, elem.second.first, elem.second.second));
    }
}

void MakeOverlapKnots(fix2knot_t &hashed_merge_knots, const contigVect_t &ctg_vect, const bool stranded, const size_t i_ovlp)
{
    for (size_t i_ctg(0); i_ctg < ctg_vect.size(); ++i_ctg)
    {
        const std::string &seq = ctg_vect[i_ctg]->GetSeq();
        uint64_t prefix = Seq2Int(seq, i_ovlp, true), suffix = Seq2Int(seq.substr(seq.size() - i_ovlp), i_ovlp, true);
        bool is_prefix_rc(false), is_suffix_rc(false);
        if (!stranded)
        {
            uint64_t prefix_rc = GetRC(prefix, i_ovlp), suffix_rc = GetRC(suffix, i_ovlp);
            if (prefix_rc < prefix)
            {
                prefix = prefix_rc;
                is_prefix_rc = true;
            }
            if (suffix_rc < suffix)
            {
                suffix = suffix_rc;
                is_suffix_rc = true;
            }
        }
        hashed_merge_knots.insert({prefix, MergeKnot()}).first->second.AddContig(i_ctg, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
        hashed_merge_knots.insert({suffix, MergeKnot()}).first->second.AddContig(i_ctg, is_suffix_rc, (is_suffix_rc ? "succ" : "pred"));
    }
}

const bool IsFirstContigRep(const float ctg_val1, const float ctg_val2, const std::string &rep_mode)
{
    if (rep_mode == "min")
    {
        return (ctg_val1 <= ctg_val2);
    }
    else if (rep_mode == "max")
    {
        return (ctg_val1 >= ctg_val2);
    }
    else if (rep_mode == "minabs")
    {
        return (fabs(ctg_val1) <= fabs(ctg_val2));
    }
    else // maxabs
    {
        return (fabs(ctg_val1) >= fabs(ctg_val2));
    }
}

const bool DoExtension(contigVect_t &ctg_vect, const fix2knot_t &hashed_mergeknot_list, const size_t i_ovlp,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_mat, const size_t nb_smp, const std::string &rep_mode)
{
    static std::vector<float> pred_counts, succ_counts;
    bool has_new_extensions(false);
    for (auto it = hashed_mergeknot_list.cbegin(); it != hashed_mergeknot_list.cend(); ++it)
    {
        if (!it->second.IsMergeable())
        {
            continue;
        }
        auto &pred_ctg = ctg_vect[it->second.GetSerial("pred")], &succ_ctg = ctg_vect[it->second.GetSerial("succ")];
        if (pred_ctg == nullptr || succ_ctg == nullptr)
        {
            continue;
        }
        const bool pred_rc = it->second.IsRC("pred"), succ_rc = it->second.IsRC("succ");
        if (interv_method == "pearson" &&
            CalcPearsonDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                            GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "spearman" &&
                 CalcSpearmanDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                                  GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "mac" &&
                 CalcMACDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                             GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        // the base contig should have minimum p-value or input order //
        if (IsFirstContigRep(pred_ctg->GetRepVal(), succ_ctg->GetRepVal(), rep_mode)) // merge right to left
        {
            if (pred_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                pred_ctg->LeftExtend(std::move(succ_ctg), !succ_rc, i_ovlp); // reverse this ctg and extend at right <=> reverse right ctg and extend to left of this ctg
            }
            else
            {
                pred_ctg->RightExtend(std::move(succ_ctg), succ_rc, i_ovlp);
            }
        }
        else // merge left to right
        {
            if (succ_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                succ_ctg->RightExtend(std::move(pred_ctg), !pred_rc, i_ovlp); // reverse this ctg and extend at left <=> reverse left ctg and extend to right of this ctg
            }
            else
            {
                succ_ctg->LeftExtend(std::move(pred_ctg), pred_rc, i_ovlp);
            }
        }
        has_new_extensions = true;
    }
    return has_new_extensions;
}

void Int2Seq(std::string &seq, const uint64_t code, const size_t k_length);
void PrintMergeKnots(const fix2knot_t &hashed_merge_knots, const contigVect_t &ctg_vect, const size_t k_len)
{
    for (const auto &elem : hashed_merge_knots)
    {
        if (elem.second.IsMergeable())
        {
            std::string fix,
                contig_pred = ctg_vect[elem.second.GetSerial("pred")]->GetSeq(),
                contig_succ = ctg_vect[elem.second.GetSerial("succ")]->GetSeq();
            Int2Seq(fix, elem.first, k_len);
            std::cout << fix << ": " << contig_pred << " ======= " << contig_succ << std::endl;
        }
    }
}

void PrintResults(const std::vector<std::string> &colname_vect, const contigVect_t &ctg_vect,
                  const std::map<uint64_t, std::pair<size_t, float>> &code_posval_map,
                  std::ifstream &idx_mat, const size_t nb_smp, const bool has_value, const std::string &out_mode,
                  const size_t k_len, const bool stranded)
{
    if (!out_mode.empty())
    {
        std::cout << "contig";
        if (has_value)
        {
            std::cout << "\trep-value";
        }
        for (const auto &s : colname_vect)
        {
            std::cout << "\t" << s;
        }
        std::cout << std::endl;
    }
    std::vector<float> count_vect;
    std::string rep_seq;
    for (const auto &elem : ctg_vect)
    {
        std::cout << elem->GetSeq();

        if (has_value)
        {
            std::cout << "\t" << elem->GetRepVal();
        }
        std::cout << "\t" << GetTagSeq(rep_seq, idx_mat, elem->GetRepPos(), nb_smp);

        if (!out_mode.empty())
        {
            if (out_mode == "rep")
            {
                GetCountVect(count_vect, idx_mat, elem->GetRepPos(), nb_smp); // output sample count vector of representative k-mer
            }
            else if (out_mode == "mean")
            {
                GetMeanCountVect(count_vect, idx_mat, nb_smp, elem->GetSeq(), code_posval_map, k_len, stranded); // output mean sample count vector
            }
            else
            {
                GetMedianCountVect(count_vect, idx_mat, nb_smp, elem->GetSeq(), code_posval_map, k_len, stranded); // output median sample count vector
            }
            for (const float x : count_vect)
            {
                std::cout << "\t" << x;
            }
            count_vect.clear();
        }
        std::cout << std::endl;
    }
}

int MergeMain(int argc, char **argv)
{
    MergeWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, sel_path, rep_mode("min"), itv_mthd("spearman"), out_path, out_mode;
    float itv_thres(0.25);
    size_t max_ovlp(0), min_ovlp(0), nb_smp, k_len;
    bool stranded(false), has_value(false);
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, max_ovlp, min_ovlp, sel_path, rep_mode, itv_mthd, itv_thres, out_path, out_mode);
    std::vector<double> _smp_sum_vect; // _smp_sum_vect not needed in KaMRaT-merge
    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, _smp_sum_vect, idx_dir + "/idx-meta.bin");
    if (k_len == 0)
    {
        throw std::domain_error("KaMRaT-merge relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }
    if (k_len <= max_ovlp)
    {
        throw std::invalid_argument("max overlap (" + std::to_string(max_ovlp) + ") should not exceed k-mer length (" + std::to_string(k_len) + ")");
    }
    PrintRunInfo(idx_dir, k_len, max_ovlp, min_ovlp, stranded, sel_path, rep_mode, itv_mthd, itv_thres, out_path, out_mode);
    if (out_mode == "mean")
    {
        std::cerr << BOLDYELLOW << "[warning]" << RESET << " estimate mean counts of contigs may introduce bias" << std::endl
                  << std::endl;
    }

    std::unordered_map<uint64_t, float> sel_code_val_map;         // k-mer code and value from given list
    std::map<uint64_t, std::pair<size_t, float>> code_posval_map; // k-mer code and position from index
    if (!sel_path.empty())
    {
        has_value = LoadSelect(sel_code_val_map, sel_path, rep_mode, k_len, stranded);
    }
    LoadCodePosValMap(code_posval_map, sel_code_val_map, idx_dir + "/idx-pos.bin");
    // for (const auto &elem : code_posval_map)
    // {
    //     std::cout << elem.first << "\t" << elem.second.first << "\t" << elem.second.second << std::endl;
    // }

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

    contigVect_t ctg_vect;
    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    InitializeContigList(ctg_vect, idx_mat, nb_smp, code_posval_map);

    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    fix2knot_t hashed_merge_knots;
    for (size_t i_ovlp(max_ovlp); i_ovlp >= min_ovlp; --i_ovlp)
    {
        std::cerr << "Merging contigs with overlap " << i_ovlp << std::endl;
        bool has_new_extensions(true);
        while (has_new_extensions)
        {
            std::cerr << "\tcontig list size: " << ctg_vect.size() << std::endl;
            MakeOverlapKnots(hashed_merge_knots, ctg_vect, stranded, i_ovlp);
            // PrintMergeKnots(hashed_merge_knots, ctg_vect, k_len);
            has_new_extensions = DoExtension(ctg_vect, hashed_merge_knots, i_ovlp, itv_mthd, itv_thres, idx_mat, nb_smp, rep_mode);
            // PrintContigList(ctg_vect, kmer_count_tab, k_len, stranded, min_ovl, interv_method, interv_thres, quant_mode, idx_mat);
            ctg_vect.erase(std::remove_if(ctg_vect.begin(), ctg_vect.end(), [](const auto &elem) { return elem == nullptr; }),
                           ctg_vect.end());
            hashed_merge_knots.clear();
        }
    }
    PrintResults(colname_vect, ctg_vect, code_posval_map, idx_mat, nb_smp, has_value, out_mode, k_len, stranded);
    idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}